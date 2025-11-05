"""
Utility methods used across the workflows
"""
import hashlib
from hailtop.batch.job import BashJob
from enum import Enum
from cpg_flow import targets
from cpg_flow.utils import logger
from cpg_utils import Path
from cpg_utils.cloud import read_secret
from cpg_utils.config import ConfigError, config_retrieve, reference_path, image_path
from os.path import join
from random import randint
from typing import Any

from functools import cache

from hailtop.batch.job import Job

from cpg_utils.cromwell import CromwellOutputType, run_cromwell_workflow_from_repo_and_get_outputs
from cpg_utils.hail_batch import command, get_batch

GATK_SV_COMMIT = config_retrieve(['workflow', 'gatk_sv_commit'])


class CromwellJobSizes(Enum):
    """
    Enum for polling intervals
    """

    SMALL = 'small'
    MEDIUM = 'medium'
    LARGE = 'large'


def get_dataset_name(dataset: str) -> str:
    """
    Add -test suffix to dataset name if in test mode.
    """
    test = config_retrieve(['workflow', 'access_level']) == 'test'
    return dataset + '-test' if test else dataset


def get_dataset_names(datasets: str | list[str]) -> list[str]:
    """
    Add -test suffix to dataset names if in test mode.
    """
    return [get_dataset_name(dataset) for dataset in datasets]


def get_query_filter_from_config(field_name: str, make_tuple = True) -> tuple[str] | list[str] | None:
    """
    Get values for the specified field from the workflow.lrs_annotation config dictionary.

    Returns tuples by default, because they are needed for use with cached functions.
    """
    if values := config_retrieve(
        ['workflow', 'query_filter', field_name],
        default=None,
    ):
        if isinstance(values, str):
            values = (values,)
        if not make_tuple:
            return list(values)
        return tuple(values)
    return None


def get_sg_hash(sequencing_group_ids: list[str]) -> str:
    """
    Unique hash string of sequencing group IDs
    Reimplemented from CPG flow core for use with LRS annotation pipeline, which sometimes
    uses subsets of sequencing groups in a given multicohort.
    """
    s = ' '.join(sorted(sequencing_group_ids))
    # use a short hash to avoid exceeding the 38 character limit for Hail Batch job
    h = hashlib.sha256(s.encode()).hexdigest()[:38]
    return f'{h}_{len(sequencing_group_ids)}'


@cache
def get_family_sequencing_groups(dataset: targets.Dataset) -> dict | None:
    """
    Get the subset of sequencing groups that are in the specified families for a dataset
    Returns a dict containing the sequencing groups and a name suffix for the outputs
    """
    if not config_retrieve(['workflow', dataset.name, 'only_families'], []):
        return None
    only_family_ids = set(config_retrieve(['workflow', dataset.name, 'only_families'], []))
    # keep only the SG IDs for the families in the only_families list
    logger.info(f'Finding sequencing groups for families {only_family_ids} in dataset {dataset.name}')
    family_sg_ids = [sg.id for sg in dataset.get_sequencing_groups() if sg.pedigree.fam_id in only_family_ids]
    if not family_sg_ids:
        raise ValueError(f'No sequencing groups found for families {only_family_ids} in dataset {dataset.name}.')
    logger.info(f'Keeping only {len(family_sg_ids)} SGs from families {len(only_family_ids)} in {dataset}:')
    logger.info(only_family_ids)
    logger.info(family_sg_ids)

    h = hashlib.sha256(''.join(sorted(family_sg_ids)).encode()).hexdigest()[:4]
    name_suffix = f'{len(family_sg_ids)}_sgs-{len(only_family_ids)}_families-{h}'

    return {'family_sg_ids': family_sg_ids, 'name_suffix': name_suffix}


def get_resource_overrides_for_job(job: BashJob, job_key: str) -> BashJob:
    """
    Get the resource overrides for a job from the workflow.resource_overrides config dictionary.
     e.g. {'storage_gib': 10, 'cpu_cores': 4}
    If no overrides are found, the job is returned unchanged.
    """
    def convert_to_gib(value: str | int) -> str:
        """
        Convert a value to a string with 'Gi' suffix for Gibibytes, or 2^30 bytes.
        """
        if isinstance(value, str) and value.endswith(('G', 'Gi')):
            return value
        if isinstance(value, int):
            return f'{value}Gi'
        return f'{value}Gi'

    overrides = config_retrieve(['workflow', 'resource_overrides', job_key], {})
    if not isinstance(overrides, dict):
        raise ValueError(f'Expected a dictionary for resource overrides for {job_key}, got {overrides}')

    # Disk
    if 'storage_gib' in overrides:
        job.storage(convert_to_gib(overrides['storage_gib']))

    # Memory
    if overrides.get('memory') in ('lowmem', 'standard', 'highmem'):
        job.memory(overrides['memory'])
    elif 'memory_gib' in overrides:
        job.memory(convert_to_gib(overrides['memory_gib']))

    # Cores
    if 'cpu_cores' in overrides:
        job.cpu(overrides['cpu_cores'])

    # Spot (preemptible) instance
    if 'spot' in overrides:
        job.spot(overrides['spot'])

    return job


def get_init_batch_args_for_job(job_name: str) -> dict[str, str | int]:
    """
    Finds init_batch args for a particular job from the config.
    e.g. {'worker_memory': 'highmem'} -> 'init_batch(worker_memory="highmem")' (via kwargs)
    """
    init_batch_args: dict[str, str | int] = {}
    workflow_config = config_retrieve(['workflow', job_name], {})

    # Memory parameters
    for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
        if workflow_config.get(config_key):
            init_batch_args[batch_key] = 'highmem'
    # Cores parameter
    for key in ['driver_cores', 'worker_cores']:
        if workflow_config.get(key):
            init_batch_args[key] = workflow_config[key]

    return init_batch_args


def get_intervals_from_bed(intervals_path: Path) -> list[str]:
    """
    Read genomic intervals from a bed file.
    Increment the start position of each interval by 1 to match the 1-based
    coordinate system used by GATK.

    Returns a list of interval strings in the format 'chrN:start-end'.
    """
    with intervals_path.open('r') as f:
        intervals = []
        for line in f:
            chrom, start, end = line.strip().split('\t')
            intervals.append(f'{chrom}:{int(start)+1}-{end}')
    return intervals


def joint_calling_scatter_count(sequencing_group_count: int) -> int:
    """
    Number of partitions for joint-calling jobs (GenotypeGVCFs, VQSR, VEP),
    as a function of the sequencing group number.
    """
    if scatter_count := config_retrieve(['workflow', 'scatter_count'], default=None):
        return scatter_count

    # Estimating this is challenging because GenotypeGVCFs does not scale
    # linearly with the number of genomes.
    # Values are adjusted based on experience with the actual number of genomes.
    # e.g. 1000 scatter count was too low for 3800 genomes.
    for threshold, scatter_count in {
        4000: 1400,
        3500: 1200,
        3000: 1000,
        2000: 600,
        1000: 400,
        500: 200,
        250: 100,
    }.items():
        if sequencing_group_count >= threshold:
            return scatter_count
    return 50


def write_mapping_to_file(mapping: dict[str, str], output_path: Path) -> None:
    """
    Write a mapping to a file
    """
    with output_path.open('w') as f:
        for k, v in mapping.items():
            f.write(f'{k} {v}\n')


def es_password() -> str:
    """
    Get Elasticsearch password. Moved into a separate method to simplify
    mocking in tests.
    """
    return read_secret(
        project_id=config_retrieve(['elasticsearch', 'password_project_id']),
        secret_name=config_retrieve(['elasticsearch', 'password_secret_id']),
        fail_gracefully=False,
    )

def make_job_name(
    name: str,
    sequencing_group: str | None = None,
    participant_id: str | None = None,
    dataset: str | None = None,
    part: str | None = None,
) -> str:
    """
    Extend the descriptive job name to reflect job attributes.
    """
    if sequencing_group and participant_id:
        sequencing_group = f'{sequencing_group}/{participant_id}'
    if sequencing_group and dataset:
        name = f'{dataset}/{sequencing_group}: {name}'
    elif dataset:
        name = f'{dataset}: {name}'
    if part:
        name += f', {part}'
    return name


def get_references(keys: list[str | dict[str, str]]) -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with reference file paths.
    """
    res: dict[str, str | list[str]] = {}
    for key in keys:
        # Keys can be maps (e.g. {'MakeCohortVcf.cytobands': 'cytoband'})
        if isinstance(key, dict):
            key, ref_d_key = next(iter(key.items()))  # noqa: PLW2901
        else:
            ref_d_key = key
        # e.g. GATKSVPipelineBatch.rmsk -> rmsk
        ref_d_key = ref_d_key.split('.')[-1]
        try:
            res[key] = reference_path(f'gatk_sv/{ref_d_key}')
        except KeyError:
            res[key] = reference_path(f'broad/{ref_d_key}')
        except ConfigError:
            res[key] = reference_path(f'broad/{ref_d_key}')

    return res

def get_images(keys: list[str], allow_missing=False) -> dict[str, str]:
    """
    Dict of WDL inputs with docker image paths.

    Args:
        keys (list): all the images to get
        allow_missing (bool): if False, require all query keys to be found

    Returns:
        dict of image keys to image paths
        or AssertionError
    """
    image_keys = config_retrieve(['images']).keys()

    if not allow_missing:
        query_keys = set(keys)
        if not query_keys.issubset(image_keys):
            raise KeyError(f'Unknown image keys: {query_keys - image_keys}')

    return {k: image_path(k) for k in image_keys if k in keys}

@cache
def create_polling_intervals() -> dict:
    """
    Set polling intervals for cromwell status
    these values are integers, indicating seconds
    for each job size, there is a min and max value
    analysis-runner implements a backoff-retrier when checking for
    success, with a minimum value, gradually reaching a max ceiling

    a config section containing overrides would look like

    [cromwell_polling_intervals.medium]
    min = 69
    max = 420
    """

    # create this dict with default values
    polling_interval_dict = {
        CromwellJobSizes.SMALL: {'min': 30, 'max': 140},
        CromwellJobSizes.MEDIUM: {'min': 40, 'max': 400},
        CromwellJobSizes.LARGE: {'min': 200, 'max': 2000},
    }

    # update if these exist in config
    for job_size in CromwellJobSizes:
        if val := config_retrieve(['cromwell_polling_intervals', job_size.value], False):
            polling_interval_dict[job_size].update(val)
    return polling_interval_dict


def add_gatk_sv_jobs(
    dataset_name: str,
    wfl_name: str,
    # "dict" is invariant (supports updating), "Mapping" is covariant (read-only)
    # we have to support inputs of type dict[str, str], so using Mapping here:
    input_dict: dict[str, Any],
    expected_out_dict: dict[str, Path | list[Path]],
    sequencing_group_id: str | None = None,
    driver_image: str | None = None,
    labels: dict[str, str] | None = None,
    job_size: CromwellJobSizes = CromwellJobSizes.MEDIUM,
) -> list[Job]:
    """
    Generic function to add a job that would run one GATK-SV workflow.
    """

    # create/retrieve dictionary of polling intervals for cromwell status
    polling_intervals = create_polling_intervals()

    # obtain upper and lower polling bounds for this job size
    polling_minimum = randint(polling_intervals[job_size]['min'], polling_intervals[job_size]['min'] * 2)  # noqa: S311
    polling_maximum = randint(polling_intervals[job_size]['max'], polling_intervals[job_size]['max'] * 2)  # noqa: S311

    # If a config section exists for this workflow, apply overrides
    if override := config_retrieve(['resource_overrides', wfl_name], False):
        input_dict |= override

    # Where Cromwell writes the output.
    # Will be different from paths in expected_out_dict:
    output_prefix = f'gatk_sv/output/{wfl_name}/{dataset_name}'
    if sequencing_group_id:
        output_prefix = join(output_prefix, sequencing_group_id)

    outputs_to_collect: dict[str, CromwellOutputType] = {}
    for key, value in expected_out_dict.items():
        if isinstance(value, list):
            outputs_to_collect[key] = CromwellOutputType.array_path(name=f'{wfl_name}.{key}', length=len(value))
        else:
            outputs_to_collect[key] = CromwellOutputType.single_path(f'{wfl_name}.{key}')

    driver_image = driver_image or image_path('cpg_workflows')

    # pre-process input_dict
    paths_as_strings: dict = {}
    for key, value in input_dict.items():
        if isinstance(value, Path):
            paths_as_strings[f'{wfl_name}.{key}'] = str(value)
        elif isinstance(value, list | set):
            paths_as_strings[f'{wfl_name}.{key}'] = [str(v) for v in value]
        else:
            paths_as_strings[f'{wfl_name}.{key}'] = value

    job_prefix = make_job_name(wfl_name, sequencing_group=sequencing_group_id, dataset=dataset_name)

    submit_j, output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
        b=get_batch(),
        job_prefix=job_prefix,
        dataset=config_retrieve(['workflow', 'dataset']),
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow=f'{wfl_name}.wdl',
        libs=['.'],
        output_prefix=output_prefix,
        input_dict=paths_as_strings,
        outputs_to_collect=outputs_to_collect,
        driver_image=driver_image,
        copy_outputs_to_gcp=config_retrieve(['workflow', 'copy_outputs'], False),
        labels=labels,
        min_watch_poll_interval=polling_minimum,
        max_watch_poll_interval=polling_maximum,
        time_limit_seconds=config_retrieve(['workflow', 'time_limit_seconds'], None),
    )

    copy_j = get_batch().new_job(f'{job_prefix}: copy outputs')
    copy_j.image(driver_image)
    cmds = []
    for key, resource in output_dict.items():
        out_path = expected_out_dict[key]
        if isinstance(resource, list):
            for source, dest in zip(resource, out_path, strict=False):
                cmds.append(f'gcloud storage cp "$(cat {source})" "{dest}"')
        else:
            cmds.append(f'gcloud storage cp "$(cat {resource})" "{out_path}"')
    copy_j.command(command(cmds, setup_gcp=True))
    return [submit_j, copy_j]
