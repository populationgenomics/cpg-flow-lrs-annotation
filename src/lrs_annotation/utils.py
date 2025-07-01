"""
suggested location for any utility methods or constants used across multiple stages
"""

from datetime import datetime
from functools import cache

from loguru import logger

from cpg_utils import Path
from cpg_utils.cloud import read_secret
from cpg_utils.config import config_retrieve
from metamist.graphql import gql, query

DATE_STRING: str = datetime.now().strftime('%y-%m-%d')  # noqa: DTZ005

VCF_QUERY = gql(
    """
    query MyQuery($dataset: String!, $seqTypes: [String!], $analysisType: String!, $metaFilter: JSON) {
      project(name: $dataset) {
        sequencingGroups(type: {in_: $seqTypes}, technology: {eq: "long-read"}) {
          id
          type
          technology
          platform
          analyses(type: {eq: $analysisType}, meta: {filter: $metaFilter}) {
            meta
            output
            outputs
          }
        }
      }
    }
    """,
)

LRS_IDS_QUERY = gql(
    """
    query MyQuery($dataset: String!, $seqTypes: [String!]) {
    project(name: $dataset) {
      sequencingGroups(technology: {eq: "long-read"}, type: {in_: $seqTypes}) {
          id
          sample {
            externalId
            meta
            participant {
              externalId
              reportedSex
            }
          }
        }
      }
    }
    """,
)


def get_config_options_as_tuple(field_name: str) -> tuple[str] | None:
    """
    Get a tuple of values from the config, for use with cached functions
    """
    if values := config_retrieve(
        ['workflow', 'lrs_annotation', field_name],
        default=None,
    ):
        if isinstance(values, str):
            values = (values,)
        return tuple(values)
    return None


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
    if scatter_count := config_retrieve(['workflow', 'scatter_count']):
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


def write_mapping_to_file(mapping: dict[str, str], output_file: Path) -> None:
    """
    Write a mapping to a file
    """
    with open(output_file, 'w') as f:
        for key, value in mapping.items():
            f.write(f'{key}\t{value}\n')


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


def find_sgs_to_skip(sg_vcf_dict: dict[str, dict]) -> set[str]:
    """
    Find the SGs to skip in the reformatting stage
    These are the parents if the family has been joint-called
    """
    sgs_to_skip = set()
    joint_called_families = set()
    for vcf_analysis in sg_vcf_dict.values():
        analysis_meta = vcf_analysis['meta']
        if analysis_meta.get('joint_called', False):
            joint_called_families.add(analysis_meta.get('family_id', ''))
    for sg_id, vcf_analysis in sg_vcf_dict.items():
        analysis_meta = vcf_analysis['meta']
        # Skip the parents if the family has been joint-called
        # Parents are identified by their participant ID ending in 01 or 02
        if (
            analysis_meta.get('family_id', '') in joint_called_families
            and analysis_meta.get('participant_id', '').endswith(('01', '02'))
        ):
            sgs_to_skip.add(sg_id)
    return sgs_to_skip


@cache
def query_for_lrs_vcfs(
    dataset_name: str,
    sequencing_types: tuple[str],
    variant_types: tuple[str],
    variant_callers: tuple[str] | None,
    pipeface_versions: tuple[str] | None,
    joint_called: bool = False,
    verbose: bool = False,
) -> dict[str, dict]:
    """
    Query metamist for the long-read sequencing VCFs
    Return a dictionary of each CPG ID and its corresponding VCF
    This result is cached - used in a SequencingGroupStage, but we only want to query for it once instead of once/SG
    Uses tuples for the inputs, so that they can be cached

    Args:
        dataset_name (str): the name of the dataset
        sequencing_types (tuple[str]): allowed sequencing types, e.g. ('genome', 'adaptive_sampling')
        variant_types (tuple[str] | None): allowed variant types, e.g. ('SNV', 'SV')
        variant_callers: (tuple[str] | None): allowed variant callers, e.g. ('deeptrio', 'sniffles')
        pipeface_versions (tuple[str] | None): allowed pipeface versions, e.g. ('v0.6.1', 'v0.7.0')
        joint_called (bool): whether to look for and preference joint-called VCFs
        verbose (bool): whether to print skipped VCFs to the log

    Returns:
        a dictionary of the SG IDs and their corresponding VCFs
    """
    single_sample_vcfs: dict[str, dict] = {}

    meta_filter = {
        'variant_type': {'in_': variant_types} if variant_types else None,
        'caller': {'in_': variant_callers} if variant_callers else None,
        'pipeface_version': {'in_': pipeface_versions} if pipeface_versions else None,
        'joint_called': False,
    }
    single_sample_vcfs_query_results = query(
        VCF_QUERY,
        variables={
            'dataset': dataset_name,
            'seqTypes': sequencing_types,
            'analysisType': 'variant_calling',
            'metaFilter': meta_filter,
        },
    )
    for sg in single_sample_vcfs_query_results['project']['sequencingGroups']:
        for analysis in sg['analyses']:
            single_sample_vcfs[sg['id']] = {
                'output': analysis['output'],
                'meta': analysis['meta'],
            }

    if not joint_called:
        return single_sample_vcfs

    # If joint_called is True, we need to query for the joint-called VCFs
    # And prefer them over the single-sample VCFs
    if verbose:
        logger.info(
            f'Finding {variant_types} joint-called VCFs in {dataset_name} for sequencing types: {sequencing_types},'
            f' callers: {variant_callers}, pipeface versions: {pipeface_versions}',
        )
    joint_called_vcfs: dict[str, dict] = {}
    meta_filter['joint_called'] = True
    joint_called_vcfs_query_results = query(
        VCF_QUERY,
        variables={
            'dataset': dataset_name,
            'seqTypes': sequencing_types,
            'analysisType': 'joint_variant_calling',
            'metaFilter': meta_filter
        },
    )
    for sg in joint_called_vcfs_query_results['project']['sequencingGroups']:
        for analysis in sg['analyses']:
            joint_called_vcfs[sg['id']] = {
                'output': analysis['output'],
                'meta': analysis['meta'],
            }

    # Prefer the joint-called VCFs over the single-sample VCFs
    sg_vcfs = {}
    for sg_id, single_sample_vcf in single_sample_vcfs.items():
        if sg_id not in joint_called_vcfs:
            sg_vcfs[sg_id] = single_sample_vcf
            continue
        sg_vcfs[sg_id] = joint_called_vcfs[sg_id]
    for sg_id, joint_called_vcf in joint_called_vcfs.items():
        if sg_id not in sg_vcfs:
            sg_vcfs[sg_id] = joint_called_vcf

    # Remove the parents entries if their family has a joint-called VCF
    sgs_to_skip = find_sgs_to_skip(sg_vcfs)
    return_dict = {}
    for sg_id, vcf_analysis in sg_vcfs.items():
        if sg_id in sgs_to_skip:
            logger.info(f'Skipping {sg_id} as it is a parent in a joint-called VCF')
            continue
        return_dict[sg_id] = vcf_analysis

    return return_dict


def query_for_lrs_to_sg_id_and_sex_mapping(datasets: list[str]):
    """
    Query metamist for the LRS ID corresponding to each sequencing group ID, and to its participant's sex
    """
    lrs_sgid_mapping = {}
    lrs_id_sex_mapping = {}
    for dataset in datasets:
        query_results = query(LRS_IDS_QUERY, variables={'dataset': dataset})
        for sg in query_results['project']['sequencingGroups']:
            sample = sg['sample']
            participant = sample['participant']
            lrs_id = sample['meta'].get('lrs_id', None)
            if not lrs_id:
                logger.warning(
                    f'{dataset} :: No LRS ID found for {participant["externalId"]} - {sample["externalId"]}',
                )
                continue
            lrs_sgid_mapping[lrs_id] = sg['id']
            lrs_id_sex_mapping[lrs_id] = participant['reportedSex']
    return lrs_sgid_mapping, lrs_id_sex_mapping
