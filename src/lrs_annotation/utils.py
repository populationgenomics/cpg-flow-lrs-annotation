"""
Utility methods used across the workflows
"""
from hailtop.batch.job import BashJob
from cpg_utils import Path
from cpg_utils.cloud import read_secret
from cpg_utils.config import config_retrieve


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
        job.storage(convert_to_gib(overrides['storage_gb']))

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


def get_init_batch_args_for_job(job_name: str) -> str:
    """
    Finds init_batch args for a particular job from the config.

    Converts the dict into a string of key value pairs so it can be passed to the batch job
    via the command line.
    e.g. {'worker_memory': 'highmem'} -> 'init_batch(worker_memory="highmem")'
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

    # translate any input arguments into an embeddable String
    return ', '.join(f'{k}={v!r}' for k, v in init_batch_args.items()) if init_batch_args else ''


def parse_init_batch_args(init_batch_args: str | None) -> dict[str, str]:
    """
    Parse the init_batch_args string into a dictionary.
    e.g. 'worker_memory="highmem", driver_memory="highmem"' -> {'worker_memory': 'highmem', 'driver_memory': 'highmem'}
    """
    if not init_batch_args:
        return {}
    args = {}
    for arg in init_batch_args.split(','):
        key, value = arg.split('=')
        args[key.strip()] = value.strip().strip("'")
    return args


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
