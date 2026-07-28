from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import Batch

from lrs_annotation.scripts import mt_to_es
from lrs_annotation.utils import get_resource_overrides_for_job


def export_mt_to_elasticsearch(
    batch: Batch,
    dataset: str,
    mt_path: str,
    index_name: str,
    flag_name: str,
    req_storage: int,
    sg_ids: list[str],
    seqr_dataset_type: str,
    input_vcfs_file_path: str,
    job_name: str,
    job_attrs: dict[str, str],
):
    """
    Export the annotated Matrixtable to ElasticSearch.
    """
    es_export_job = batch.new_job(job_name, attributes=job_attrs)
    es_export_job.image(config_retrieve(['workflow', 'driver_image']))
    es_export_job.storage(f'{req_storage}Gi')
    es_export_job.spot(is_spot=config_retrieve(['elasticsearch', 'spot_instance'], default=True))
    es_export_job = get_resource_overrides_for_job(es_export_job, 'export_mt_to_elasticsearch')

    # localise the MT
    mt_name = mt_path.split('/')[-1]
    es_export_job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')

    # run the export from the localised MT - this job writes no new data, just transforms and exports over network
    es_export_job.command(
        f"""
        python3 {mt_to_es.__file__} \\
        --dataset {dataset} \\
        --mt_path $BATCH_TMPDIR/{mt_name} \\
        --index {index_name} \\
        --flag {flag_name} \\
        --sg_ids {' '.join(sg_ids)} \\
        --seqr_dataset_type {seqr_dataset_type} \\
        --path_to_input_vcfs_file {input_vcfs_file_path}
        """
    )

    return es_export_job
