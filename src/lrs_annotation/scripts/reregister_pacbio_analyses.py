import argparse

from cpg_utils import to_path
from loguru import logger
from metamist.apis import AnalysisApi
from metamist.graphql import gql, query
from metamist.models import Analysis, AnalysisStatus, AnalysisUpdateModel

PACBIO_SGS_QUERY = gql(
    """
    query pacbioSgs($dataset: String!) {
        project(name: $dataset) {
            sequencingGroups(type: {eq_: "genome"}, technology: {eq_: "long-read"}, platform: {eq_: "pacbio"}) {
                id
                analyses(type: {in_: ["pacbio_cram", "cram", "pacbio_vcf", "vcf"]}, status: {eq_: COMPLETED}) {
                    id
                    type
                    output
                    timestampCompleted
                    meta
                }
            }
        }
    }
    """
)


def main():
    """
    Find long-read PacBio SGs with "pacbio_cram" and "pacbio_vcf" analyses and re-register them as "cram" and "vcf"
    analyses respectively. This is needed because the analysis type was changed in the workflow, but the old analyses
    are still registered in Metamist.
    """
    parser = argparse.ArgumentParser(description='Register analyses for CRAM files')
    parser.add_argument('dataset', type=str, help='Dataset name to register analyses for')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    # Query for sequencing groups and their analyses
    result = query(PACBIO_SGS_QUERY, variables={'dataset': args.dataset})
    sequencing_groups = result['project']['sequencingGroups']

    analysis_api = AnalysisApi()

    sg_analyses_template = {
        'pacbio_cram': None,
        'pacbio_vcf': None,
        'cram': None,
        'vcf': None,
    }

    # dict mapping sg_id to the template above, with the analyses filled in for the analyses that exist for that SG
    sg_analyses = {}

    for sg in sequencing_groups:
        sg_id = sg['id']
        analyses = sg['analyses']
        sg_analyses[sg_id] = sg_analyses_template.copy()
        for analysis in analyses:
            if analysis['type'] in sg_analyses[sg_id]:
                sg_analyses[sg_id][analysis['type']] = analysis

        # If the "pacbio_cram" analysis exists but the "cram" analysis does not, register the "cram" analysis
        if sg_analyses[sg_id]['pacbio_cram'] and not sg_analyses[sg_id]['cram']:
            pacbio_cram_path = to_path(sg_analyses[sg_id]['pacbio_cram']['output'])
            pacbio_cram_timestamp_completed = sg_analyses[sg_id]['pacbio_cram']['timestampCompleted']
            pacbio_cram_meta = sg_analyses[sg_id]['pacbio_cram']['meta']
            if not pacbio_cram_path.exists():
                logger.info(
                    f'CRAM file for sequencing group {sg_id} does not exist at expected path: {pacbio_cram_path}'
                )
                continue
            if not args.dry_run:
                logger.info(
                    f'Registering "cram" analysis for sequencing group {sg_id} with CRAM path: {pacbio_cram_path}'
                )
                analysis_api.create_analysis(
                    project=args.dataset,
                    analysis=Analysis(
                        type='cram',
                        status=AnalysisStatus('completed'),
                        output=str(pacbio_cram_path),
                        sequencing_group_ids=[sg_id],
                        timestamp_completed=pacbio_cram_timestamp_completed,
                        meta=pacbio_cram_meta,
                    ),
                )
                # deactivate the old "pacbio_cram" analysis so it doesn't show up in the UI
                logger.info(f'Deactivating old "pacbio_cram" analysis for sequencing group {sg_id}')
                old_pacbio_cram_analysis_id = sg_analyses[sg_id]['pacbio_cram']['id']
                analysis_api.update_analysis(
                    analysis_id=old_pacbio_cram_analysis_id,
                    analysis_update_model=AnalysisUpdateModel(
                        status=AnalysisStatus('completed'),
                        active=False,
                    ),
                )
            else:
                logger.info(f'Dry run: would register new "cram" analysis for sequencing group {sg_id}:')
                logger.info(
                    f' {pacbio_cram_path} with timestamp {pacbio_cram_timestamp_completed} and meta {pacbio_cram_meta}'
                )
                logger.info(
                    f'Dry run: would deactivate old "pacbio_cram" analysis {old_pacbio_cram_analysis_id} '
                    f'for sequencing group {sg_id}'
                )

        # If the "pacbio_vcf" analysis exists but the "vcf" analysis does not, register the "vcf" analysis
        if sg_analyses['pacbio_vcf'] and not sg_analyses['vcf']:
            pacbio_vcf_path = to_path(sg_analyses[sg_id]['pacbio_vcf']['output'])
            pacbio_vcf_timestamp_completed = sg_analyses[sg_id]['pacbio_vcf']['timestampCompleted']
            pacbio_vcf_meta = sg_analyses[sg_id]['pacbio_vcf']['meta']
            if not pacbio_vcf_path.exists():
                logger.info(f'VCF file for sequencing group {sg_id} does not exist at expected path: {pacbio_vcf_path}')
                continue
            if not args.dry_run:
                logger.info(f'Registering "vcf" analysis for sequencing group {sg_id} with VCF path: {pacbio_vcf_path}')
                analysis_api.create_analysis(
                    project=args.dataset,
                    analysis=Analysis(
                        type='vcf',
                        status=AnalysisStatus('completed'),
                        output=str(pacbio_vcf_path),
                        sequencing_group_ids=[sg_id],
                        timestamp_completed=pacbio_vcf_timestamp_completed,
                        meta=pacbio_vcf_meta,
                    ),
                )
                # deactivate the old "pacbio_vcf" analysis so it doesn't show up in the UI
                logger.info(f'Deactivating old "pacbio_vcf" analysis  for sequencing group {sg_id}')
                old_pacbio_vcf_analysis_id = sg_analyses[sg_id]['pacbio_vcf']['id']
                analysis_api.update_analysis(
                    analysis_id=old_pacbio_vcf_analysis_id,
                    analysis_update_model=AnalysisUpdateModel(
                        status=AnalysisStatus('completed'),
                        active=False,
                    ),
                )
            else:
                logger.info(f'Dry run: would register new "vcf" analysis for sequencing group {sg_id}:')
                logger.info(
                    f' {pacbio_vcf_path} with timestamp {pacbio_vcf_timestamp_completed} and meta {pacbio_vcf_meta}'
                )
                logger.info(
                    f'Dry run: would deactivate old "pacbio_vcf" analysis {old_pacbio_vcf_analysis_id} '
                    f'for sequencing group {sg_id}'
                )


if __name__ == '__main__':
    main()
