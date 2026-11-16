"""
Finds CRAM files in the /long_read/long_read/cram directories of the
input datasets and corrects the paths and analysis records.
"""

import argparse
import random
import time

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
            cram_path = to_path(outputs['path'])
            bucket = storage_client.bucket(cram_path.bucket)

            cram_blob = bucket.blob(cram_path.blob)
            new_cram_str = cram_path.blob.replace('long_read/long_read/cram', 'long_read/cram')

            # Replace the extension in the cram path to derive the index and somalier paths
            crai_blob = bucket.blob(cram_path.blob.replace('.cram', '.cram.crai'))
            new_crai_str = new_cram_str.replace('.cram', '.cram.crai')
            somalier_blob = bucket.blob(cram_path.blob.replace('.cram', '.cram.somalier'))
            new_somalier_str = new_cram_str.replace('.cram', '.cram.somalier')

            if not dry_run:
                loguru.logger.info(f'Moving {cram_blob.name} to {new_cram_str} in bucket {bucket.name}')
                for old_blob_name, new_blob_name in [
                    (cram_blob.name, new_cram_str),
                    (crai_blob.name, new_crai_str),
                    (somalier_blob.name, new_somalier_str),
                ]:
                    copy_blob_with_retries(
                        bucket=bucket,
                        source_blob_name=old_blob_name,
                        destination_blob_name=new_blob_name,
                    )
                    loguru.logger.info(f'Copied. Deleting original: {old_blob_name} from bucket {bucket.name}')
                    bucket.delete_blob(old_blob_name)

                loguru.logger.info(
                    f'Updating analysis record for {analysis["id"]} with new path gs://{bucket.name}/{new_cram_str}'
                )
                AnalysisApi().update_analysis(
                    analysis_id=analysis['id'],
                    analysis_update_model=AnalysisUpdateModel(
                        outputs={'basename': f'gs://{bucket.name}/{new_cram_str}'}
                    ),
                )
            else:
                loguru.logger.info(f'Dry run: would move {cram_blob.name} to {new_cram_str} in {bucket.name}')
                loguru.logger.info(
                    f'Dry run: would update analysis record for {analysis["id"]} with new path gs://{bucket.name}/{new_cram_str}'
                )


def copy_blob_with_retries(
    bucket: storage.Bucket,
    source_blob_name: str,
    destination_blob_name: str,
):
    """
    Copy a blob from source to destination bucket with retry logic and proper token handling.
    Uses the rewrite API as suggested by the Google Cloud docs, as blob.copy can be unreliable
    for large files.

    Args:
        bucket: Bucket object where the cram lives
        source_blob_name: Original blob name (e.g. '/long_read/long_read/cram/SG12345.cram')
        destination_blob_name: Destination blob name (e.g. 'long_read/cram/SG12345.cram')

    Returns:
        The destination blob object
    """
    source_blob = bucket.blob(source_blob_name)
    destination_blob = bucket.blob(destination_blob_name)

    max_retries = 5
    initial_backoff = 1  # Initial backoff time in seconds
    max_backoff = 5  # Maximum backoff time in seconds

    retry_count = 0
    rewrite_token = None

    while True:
        try:
            # Attempt the rewrite operation
            rewrite_token_before = rewrite_token  # Store for comparison

            # The rewrite method returns (rewrite_token, bytes_rewritten, total_bytes)
            rewrite_token, bytes_rewritten, total_bytes = destination_blob.rewrite(
                source=source_blob,
                token=rewrite_token,
            )

            # Log progress for large files
            if total_bytes > 0:
                progress = (bytes_rewritten / total_bytes) * 100
                loguru.logger.info(
                    f'Rewrite progress: {progress:.1f}% ({bytes_rewritten}/{total_bytes} bytes)',
                )

            # If rewrite_token is None, the operation is complete
            if rewrite_token is None:
                loguru.logger.info(
                    f'Successfully copied {source_blob_name} to {destination_blob_name}',
                )
                return destination_blob

            # If token hasn't changed but operation isn't complete, something is wrong
            if rewrite_token == rewrite_token_before:
                raise Exception("Rewrite operation stalled - token didn't change")

            # Reset retry count when making progress
            retry_count = 0

        except Exception as e:
            retry_count += 1
            if retry_count > max_retries:
                loguru.logger.error(
                    f'Failed to copy blob after {max_retries} attempts. Error: {e}',
                )
                raise

            # Calculate backoff with exponential increase and jitter
            backoff = min(initial_backoff * (2 ** (retry_count - 1)), max_backoff)
            backoff_with_jitter = backoff * (0.5 + random.random())  # noqa: S311

            loguru.logger.warning(
                f'Error during blob copy: {e}. Retry {retry_count}/{max_retries} after {backoff_with_jitter:.1f}s',
            )
            time.sleep(backoff_with_jitter)


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
