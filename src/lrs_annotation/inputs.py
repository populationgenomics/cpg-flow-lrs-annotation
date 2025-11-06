"""
Methods for querying Metamist for long-read sequencing VCFs and related metadata.
"""
from functools import cache

from loguru import logger

from cpg_utils.config import config_retrieve

from metamist.graphql import gql, query

from utils import get_dataset_name

VCF_QUERY = gql(
    """
    query MyQuery(
        $dataset: String!, $seqTypes: [String!], $platforms: [String!], $analysisTypes: [String!], $metaFilter: JSON
    ) {
      project(name: $dataset) {
        sequencingGroups(type: {in_: $seqTypes}, technology: {eq: "long-read"}, platform: {in_: $platforms}) {
          id
          type
          technology
          platform
          analyses(type: {in_: $analysisTypes}, meta: $metaFilter) {
            id
            type
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
    query MyQuery($dataset: String!, $seqTypes: [String!], $platforms: [String!]) {
    project(name: $dataset) {
      sequencingGroups(technology: {eq: "long-read"}, type: {in_: $seqTypes}, platform: {in_: $platforms}) {
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


def find_sgs_to_skip(sg_vcf_dict: dict[str, dict]) -> set[str]:
    """
    Find sequencing groups to skip based on joint-called VCFs.

    Since the stages in the workflows are run per sequencing group,
    we need to skip the parents of joint-called families, even if the
    parents have their own single-sample VCFs.
    """
    sgs_to_skip = set()
    joint_called_families = set()

    # Find all families that have joint-called VCFs
    for vcf_analysis in sg_vcf_dict.values():
        if vcf_analysis['meta'].get('joint_called', False):
            if not vcf_analysis['meta'].get('family_id'):
                logger.warning(
                    f"Joint-called VCF {vcf_analysis['output']} does not have a family ID, skipping.",
                )
                continue
            # Add the family ID to the set of joint-called families
            # This will be used to skip the parents in the next loop
            joint_called_families.add(vcf_analysis['meta'].get('family_id', ''))

    # Loop through the VCFs again to find parents in joint-called families
    for sg_id, vcf_analysis in sg_vcf_dict.items():
        analysis_meta = vcf_analysis['meta']
        # Skip the parents if the family has a joint-called VCF
        # Parents are identified by their participant ID ending with '01' or '02'
        # OR some other specified pattern, e.g. 'P', 'M'
        parental_id_suffixes = tuple(
            config_retrieve(
                ['workflow', 'query_filters', 'parental_id_suffixes'],
                default=['01', '02']
            )
        )
        if (
            analysis_meta.get('family_id', '') in joint_called_families
            and analysis_meta.get('participant_id', '').endswith(parental_id_suffixes)
        ):
            sgs_to_skip.add(sg_id)
    return sgs_to_skip


@cache
def query_for_lrs_vcfs(
    dataset_name: str
) -> dict[str, dict | list[str]]:
    """
    Query metamist for the long-read sequencing VCFs. The VCFs are filtered by the values specified in
    the workflow.query_filters config dictionary, using the analysis.meta field.

    Mandatory query filters:
      - sequencing_types   - e.g [ 'genome', 'adaptive_sampling', ]
      - sequencing_platforms - e.g [ 'pacbio', 'oxford-nanopore', ]
      - analysis_types     - e.g [ 'vcf', 'pacbio_vcf', ]
      - variant_types      - e.g [ 'SNV, 'snp_indel', ]
      - variant_callers    - e.g [ 'deeptrio', 'deepvariant', ]
      - pipeface_versions  - e.g [ 'v0.8.0', ]

    If prefer_joint_called is set to True, the query will return joint-called VCFs

    This result cached since we only want to query for it once, instead of once per SG stage.

    Args:
        dataset_name (str): the name of the dataset

    Returns:
        'vcfs': dict[str, dict]
            A dictionary of sequencing group IDs and their corresponding VCFs and metadata.
        'sg_ids': list[str]
            A list of sequencing group IDs that are present in the workflow run.
        Note: Not every SG will have a VCF, since some may be parents in joint-called families.
        These SGs will still be included in the 'sg_ids' list, but their VCFs will not be present in the 'vcfs' dict.
    """
    if config_retrieve(['workflow', 'access_level']) == 'test' and not dataset_name.endswith('-test'):
        dataset_name += '-test'
    query_filters: dict = config_retrieve(['workflow', 'query_filters'], default={})
    if not query_filters:
        raise ValueError('No query filters found in the config file. Please check your configuration.')
    for key in ['sequencing_types', 'analysis_types', 'variant_types', 'variant_callers', 'pipeface_versions']:
        if key not in query_filters:
            raise ValueError(f'Missing required query filter: {key}')

    sequencing_types = tuple(query_filters.get('sequencing_types', []))
    sequencing_platforms = tuple(query_filters.get('sequencing_platforms', []))
    analysis_types = tuple(query_filters.get('analysis_types', []))
    variant_types = tuple(query_filters.get('variant_types', []))
    variant_callers = tuple(query_filters.get('variant_callers', []))
    pipeface_versions = tuple(query_filters.get('pipeface_versions', []))

    verbose = config_retrieve(['workflow', 'verbose'], default=False)
    if verbose:
        logger.info(
            f'{dataset_name} :: Finding {variant_types} single-sample VCFs for '
            f' sequencing types: {sequencing_types}, sequence platforms: {sequencing_platforms},'
            f' analysis types: {analysis_types}, callers: {variant_callers}, pipeface versions: {pipeface_versions}',
        )
    single_sample_vcfs: dict[str, dict] = {}

    meta_filter = {
        'variant_type': {'in_': variant_types} if variant_types else None,
        'caller': {'in_': variant_callers} if variant_callers else None,
        'pipeface_version': {'in_': pipeface_versions} if pipeface_versions else None,
        'joint_called': {'eq': 'false'},
    }
    single_sample_vcfs_query_results = query(
        VCF_QUERY,
        variables={
            'dataset': dataset_name,
            'seqTypes': sequencing_types,
            'platforms': sequencing_platforms,
            'analysisTypes': analysis_types,
            'metaFilter': meta_filter,
        },
    )
    for sg in single_sample_vcfs_query_results['project']['sequencingGroups']:
        for analysis in sg['analyses']:
            single_sample_vcfs[sg['id']] = {
                'vcf': analysis['output'],
                'meta': analysis['meta'],
            }
    if verbose:
        logger.info(
            f'{dataset_name} :: Found {len(single_sample_vcfs)} single-sample VCFs.',
        )
    if not query_filters.get('prefer_joint_called', False):
        # If we are not preferring joint-called VCFs, just return the single-sample VCFs
        return {
            'sg_ids': list(single_sample_vcfs.keys()),
            'vcfs': single_sample_vcfs,
        }

    # If joint_called is True, we need to query for the joint-called trio VCFs
    # and prefer them over the single-sample VCFs of parents in joint-called trios
    if verbose:
        logger.info(
            f'Finding {variant_types} joint-called VCFs in {dataset_name} for '
            f' sequencing types: {sequencing_types}, sequence platforms: {sequencing_platforms},'
            f' analysis types: {analysis_types}, callers: {variant_callers}, pipeface versions: {pipeface_versions}',
        )
    joint_called_vcfs: dict[str, dict] = {}
    meta_filter['joint_called'] = {'eq': 'true'}
    joint_called_vcfs_query_results = query(
        VCF_QUERY,
        variables={
            'dataset': dataset_name,
            'seqTypes': sequencing_types,
            'platforms': sequencing_platforms,
            'analysisType': analysis_types,
            'metaFilter': meta_filter
        },
    )
    for sg in joint_called_vcfs_query_results['project']['sequencingGroups']:
        for analysis in sg['analyses']:
            joint_called_vcfs[sg['id']] = {
                'vcf': analysis['output'],
                'meta': analysis['meta'],
            }

    # Prefer the joint-called VCFs over the single-sample VCFs of parents in joint-called trios
    sg_vcfs = {}
    for sg_id, single_sample_vcf in single_sample_vcfs.items():
        if sg_id not in joint_called_vcfs:
            sg_vcfs[sg_id] = single_sample_vcf
            continue
        sg_vcfs[sg_id] = joint_called_vcfs[sg_id]
    for sg_id, joint_called_vcf in joint_called_vcfs.items():
        if sg_id not in sg_vcfs:
            sg_vcfs[sg_id] = joint_called_vcf

    # Remove the parents entries if their family has a joint-called trio VCF
    sgs_to_skip = find_sgs_to_skip(sg_vcfs)
    vcfs_for_sgs = {}
    for sg_id, vcf_analysis in sg_vcfs.items():
        if sg_id in sgs_to_skip:
            if verbose:
                logger.info(f'Skipping {sg_id} as it is a parent in a joint-called VCF')
            continue
        vcfs_for_sgs[sg_id] = vcf_analysis

    return {
        # SG IDs include the skipped parental ones, so that we can still
        # reference them in the workflow, but the vcfs dict only contains
        # the VCFs that are not skipped.
        'sg_ids': list(vcfs_for_sgs.keys()) + list(sgs_to_skip),
        'vcfs': vcfs_for_sgs,
    }


def query_for_lrs_mappings(
    dataset_names: list[str],
    sequencing_types: list[str],
    sequencing_platforms: list[str],
    ) -> dict[str, dict[str, str]]:
    """
    Query metamist for the LRS ID to SG ID mapping, and it's associated participant's sex
    """
    lrs_mappings = {}
    for dataset in dataset_names:
        query_results = query(
            LRS_IDS_QUERY,
            variables={'dataset': dataset, 'seqTypes': sequencing_types, 'platforms': sequencing_platforms}
        )
        for sg in query_results['project']['sequencingGroups']:
            sample = sg['sample']
            participant = sample['participant']
            lrs_id = sample['meta'].get('lrs_id', None)
            if not lrs_id:
                logger.warning(
                    f'{dataset} :: No LRS ID found for {participant["externalId"]} - {sample["externalId"]}',
                )
                continue
            lrs_mappings[lrs_id] = {'sg_id': sg['id'], 'sex': participant['reportedSex']}
    return lrs_mappings


def get_sgs_from_datasets(multicohort_datasets: list[str]) -> dict[str, list[str] | dict]:
    """
    Returns the sequencing group IDs from multicohort datasets, filtered to the
    sequencing groups that are actually present in the workflow run.
    """
    sg_ids: list[str] = []
    vcfs: dict[str, dict] = {}
    for dataset in multicohort_datasets:
        sg_ids.extend(query_for_lrs_vcfs(get_dataset_name(dataset))['sg_ids'])
        vcfs.update(query_for_lrs_vcfs(get_dataset_name(dataset))['vcfs'])  # type: ignore[arg-type]
    return {
        'sg_ids': sg_ids,
        'vcfs': vcfs
    }
