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
            sequencingGroups(type: {eq: "genome"}, technology: {eq: "long-read"}, platform: {eq: "pacbio"}) {
                id
                analyses(type: {in_: ["pacbio_cram", "cram", "pacbio_vcf", "vcf"]}, status: {eq: COMPLETED}) {
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


def main():  # noqa: PLR0915
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
        'pacbio_cram': [],
        'pacbio_vcf': [],
        'cram': [],
        'vcf': [],
    }

    # dict mapping sg_id to the template above, with the analyses filled in for the analyses that exist for that SG
    sg_analyses = {}

    for sg in sequencing_groups:
        sg_id = sg['id']
        analyses = sg['analyses']
        sg_analyses[sg_id] = sg_analyses_template.copy()
        for analysis in analyses:
            if analysis['type'] in sg_analyses[sg_id]:
                if not sg_analyses[sg_id][analysis['type']]:
                    sg_analyses[sg_id][analysis['type']] = []
                sg_analyses[sg_id][analysis['type']].append(analysis)

        # expect one cram / pacbio_cram, but multiple vcf / pacbio_vcf analyses (e.g. from different variant callers)
        if len(sg_analyses[sg_id]['pacbio_cram']) == 1:
            sg_analyses[sg_id]['pacbio_cram'] = sg_analyses[sg_id]['pacbio_cram'][0]
        elif len(sg_analyses[sg_id]['pacbio_cram']) > 1:
            logger.warning(f'SG {sg_id} has multiple "pacbio_cram" analyses, which is unexpected. Skipping this SG.')

        if len(sg_analyses[sg_id]['cram']) == 1:
            sg_analyses[sg_id]['cram'] = sg_analyses[sg_id]['cram'][0]
        elif len(sg_analyses[sg_id]['cram']) > 1:
            logger.warning(f'SG {sg_id} has multiple "cram" analyses, which is unexpected. Skipping this SG.')

        # If the "pacbio_cram" analysis exists but the "cram" analysis does not, register the "cram" analysis
        # and always deactivate the old "pacbio_cram" analysis
        if sg_analyses[sg_id]['pacbio_cram']:
            pacbio_cram_path = to_path(sg_analyses[sg_id]['pacbio_cram']['output'])
            if not pacbio_cram_path.exists():
                logger.info(
                    f'CRAM file for sequencing group {sg_id} does not exist at expected path: {pacbio_cram_path}'
                )
                # try the standard location
                pacbio_cram_path = to_path(f'gs://cpg-{args.dataset}-main/long_read/cram/{sg_id}.cram')
                if not pacbio_cram_path.exists():
                    logger.info(
                        f'CRAM file for sequencing group {sg_id} does not exist at expected path: {pacbio_cram_path}'
                    )
                    logger.info(f'Skipping SG {sg_id} since the CRAM file does not exist')
                    continue
            pacbio_cram_timestamp_completed = sg_analyses[sg_id]['pacbio_cram']['timestampCompleted']
            pacbio_cram_meta = sg_analyses[sg_id]['pacbio_cram']['meta']
            pacbio_cram_analysis_id = sg_analyses[sg_id]['pacbio_cram']['id']
            if not args.dry_run:
                # If a "cram" analysis does not already exist, register it
                if not sg_analyses[sg_id]['cram']:
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
                analysis_api.update_analysis(
                    analysis_id=pacbio_cram_analysis_id,
                    analysis_update_model=AnalysisUpdateModel(
                        status=AnalysisStatus('completed'),
                        active=False,
                    ),
                )
            else:
                if not sg_analyses[sg_id]['cram']:
                    logger.info(f'Dry run: would register new "cram" analysis for sequencing group {sg_id}:')
                    logger.info(
                        f' {pacbio_cram_path} with timestamp {pacbio_cram_timestamp_completed} '
                        f'and meta {pacbio_cram_meta}'
                    )
                logger.info(
                    f'Dry run: would deactivate old "pacbio_cram" analysis {pacbio_cram_analysis_id} '
                    f'for sequencing group {sg_id}'
                )

        # If the "pacbio_vcf" analyses exist but the corresponding "vcf" analyses do not, register the "vcf" analyses
        # and always deactivate the old "pacbio_vcf" analyses
        if sg_analyses[sg_id]['pacbio_vcf']:
            for pacbio_vcf_analysis in sg_analyses[sg_id]['pacbio_vcf']:
                pacbio_vcf_path = to_path(pacbio_vcf_analysis['output'])
                pacbio_vcf_timestamp_completed = pacbio_vcf_analysis['timestampCompleted']
                pacbio_vcf_meta = pacbio_vcf_analysis['meta']
                pacbio_vcf_analysis_id = pacbio_vcf_analysis['id']
                if not pacbio_vcf_path.exists():
                    logger.info(
                        f'VCF file for sequencing group {sg_id} does not exist at expected path: {pacbio_vcf_path}'
                    )
                    continue
                if not args.dry_run:
                    # Too hard to check for perfect matches of existing "vcf" analyses so just register them regardless
                    logger.info(
                        f'Registering "vcf" analysis for sequencing group {sg_id} with VCF path: {pacbio_vcf_path}'
                    )
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
                    logger.info(f'Deactivating old "pacbio_vcf" analysis for sequencing group {sg_id}')
                    analysis_api.update_analysis(
                        analysis_id=pacbio_vcf_analysis_id,
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
                        f'Dry run: would deactivate old "pacbio_vcf" analysis {pacbio_vcf_analysis_id} '
                        f'for sequencing group {sg_id}'
                    )


if __name__ == '__main__':
    main()
