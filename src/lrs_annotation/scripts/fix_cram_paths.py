"""
Finds CRAM files in the /long_read/long_read/cram directories of the
input datasets and corrects the paths and analysis records.
"""

import argparse

import loguru
from google.cloud import storage

from cpg_utils import to_path
from metamist.apis import AnalysisApi
from metamist.graphql import gql, query
from metamist.models import AnalysisUpdateModel

CRAM_ANALYSES_QUERY = gql("""

query cramAnalyses($dataset: String!) {
    project(name: $dataset) {
        sequencingGroups(technology: {eq: "long-read"}) {
            id
            analyses(type: {eq: "CRAM"}) {
                id
                outputs
            }
        }
    }
}
""")

storage_client = storage.Client()


def fix_cram_paths(result, dry_run=False):
    """
    Find long read cram analyses at the wrong path and move them to the correct location.
    """
    for sg in result['project']['sequencingGroups']:
        for analysis in sg['analyses']:
            if 'long_read/long_read/cram' not in analysis['outputs']['path']:
                continue

            outputs = analysis['outputs']
            source_path = to_path(outputs['path'])
            bucket = storage_client.bucket(source_path.bucket)

            source_blob = bucket.blob(source_path.blob)
            new_blob_str = source_path.blob.replace('long_read/long_read/cram', 'long_read/cram')

            crai_blob = bucket.blob(source_path.blob.replace('.cram', '.cram.crai'))
            new_crai_str = new_blob_str.replace('.cram', '.cram.crai')
            somalier_blob = bucket.blob(source_path.blob.replace('.cram', '.cram.somalier'))
            new_somalier_str = new_blob_str.replace('.cram', '.cram.somalier')

            if not dry_run:
                loguru.logger.info(f'Moving {source_blob.name} to {new_blob_str} in bucket {bucket.name}')
                bucket.copy_blob(source_blob, bucket, new_blob_str)
                bucket.copy_blob(crai_blob, bucket, new_crai_str)
                bucket.copy_blob(somalier_blob, bucket, new_somalier_str)
                loguru.logger.info(f'Copied. Deleting original: {source_path.blob} from bucket {bucket.name}')
                bucket.delete_blob(source_path.blob)
                bucket.delete_blob(crai_blob.name)
                bucket.delete_blob(somalier_blob.name)
                loguru.logger.info(
                    f'Updating analysis record for {analysis["id"]} with new path gs://{bucket.name}/{new_blob_str}'
                )
                AnalysisApi().update_analysis(
                    analysis_id=analysis['id'],
                    analysis_update_model=AnalysisUpdateModel(outputs={'basename': f'gs://{bucket.name}/{new_blob_str}'}),
                )
            else:
                loguru.logger.info(f'Dry run: would move {source_blob.name} to {new_blob_str} in {bucket.name}')
                loguru.logger.info(
                    f'Dry run: would update analysis record for {analysis["id"]} with new path gs://{bucket.name}/{new_blob_str}'
                )


def main():
    parser = argparse.ArgumentParser(description='Fix CRAM paths in datasets')
    parser.add_argument('-d', '--dataset', help='Dataset(s) to process', nargs='+')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    for dataset in args.dataset:
        result = query(CRAM_ANALYSES_QUERY, {'dataset': dataset})

        fix_cram_paths(result, dry_run=args.dry_run)


if __name__ == '__main__':
    main()
