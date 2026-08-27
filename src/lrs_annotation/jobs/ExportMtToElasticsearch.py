from cpg_flow.utils import tshirt_mt_sizing
from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import Batch

from lrs_annotation.scripts import mt_to_es
from lrs_annotation.utils import get_resource_overrides_for_job


def export_mt_to_elasticsearch(
    batch: Batch,
    dataset: str,
    sg_ids: list[str],
    mt_path: str,
    index_name: str,
    done_flag: str,
    input_vcfs_file_path: Path,
    job_attrs: dict[str, str],
):
    """
    Export the annotated Matrixtable to ElasticSearch.
    """
    req_storage = tshirt_mt_sizing(
        sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
        cohort_size=len(sg_ids),
    )
    mt_name = mt_path.rsplit('/', maxsplit=1)[-1]
    job_name = f'Export {index_name} from {mt_name}'
    es_export_job = batch.new_job(job_name, attributes=job_attrs)
    es_export_job.image(config_retrieve(['workflow', 'driver_image']))
    es_export_job.storage(f'{req_storage}Gi')
    es_export_job.spot(is_spot=config_retrieve(['elasticsearch', 'spot_instance'], default=True))
    es_export_job = get_resource_overrides_for_job(es_export_job, 'export_mt_to_elasticsearch')

    # localise the MT
    es_export_job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')

    # run the export from the localised MT - this job writes no new data, just transforms and exports over network
    es_export_job.command(
        f"""
        python3 {mt_to_es.__file__} \\
        --dataset {dataset} \\
        --sg_ids {' '.join(sg_ids)} \\
        --mt_path $BATCH_TMPDIR/{mt_name} \\
        --index {index_name} \\
        --done_flag {done_flag} \\
        --path_to_input_vcfs_file {input_vcfs_file_path!s}
        """
    )

    return es_export_job
