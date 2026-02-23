from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from hailtop.batch.job import Job

from lrs_annotation.scripts.snps_indels import vcf_to_unannotated_mt


def vcf_to_unannotated_mt_job(
    vcf_path: Path,
    out_mt_path: Path,
    job_attrs: dict | None = None,
) -> Job:
    """
    Convert VCF to unannotated matrix table for Talos.
    """
    j = get_batch().new_job('Convert VCF to unannotated MT', job_attrs)
    j.image(config_retrieve(['workflow', 'driver_image']))
    j.command(
        f"""
        python3 {vcf_to_unannotated_mt.__file__} \\
            --vcf_path {vcf_path} \\
            --out_mt_path {out_mt_path}
        """
    )
    return j
