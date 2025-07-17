from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from lrs_annotation.scripts import subset_mt_to_sgs
from lrs_annotation.scripts.svs import annotate_dataset_mt


def annotate_dataset_jobs_sv(
    mt_path: Path,
    sg_ids: list[str],
    out_mt_path: Path,
    tmp_prefix: Path,
    job_attrs: dict | None = None,
    depends_on: list[Job] | None = None,
) -> list[Job]:
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr).
    """
    sample_ids_list_path = tmp_prefix / 'sample-list.txt'
    if not config_retrieve(['workflow', 'dry_run'], default=False):
        with sample_ids_list_path.open('w') as f:
            f.write(','.join(sg_ids))

    subset_mt_path = tmp_prefix / 'cohort-subset.mt'

    subset_j = get_batch().new_job('Subset cohort to dataset', (job_attrs or {}) | {'tool': 'hail query'})
    subset_j.image(config_retrieve(['workflow', 'driver_image']))
    assert sg_ids
    subset_j.command(
        f"""
        python3 {subset_mt_to_sgs.__file__} \\
            --mt_path {mt_path} \\
            --sg_ids {','.join(sg_ids)} \\
            --out_mt_path {subset_mt_path}
        """
    )
    if depends_on:
        subset_j.depends_on(*depends_on)

    annotate_j = get_batch().new_job('Annotate dataset', (job_attrs or {}) | {'tool': 'hail query'})
    annotate_j.image(config_retrieve(['workflow', 'driver_image']))
    annotate_j.command(
        f"""
        python3 {annotate_dataset_mt.__file__} \\
            --mt_path {subset_mt_path} \\
            --out_mt_path {out_mt_path}
        """
    )
    annotate_j.depends_on(subset_j)
    return [subset_j, annotate_j]
