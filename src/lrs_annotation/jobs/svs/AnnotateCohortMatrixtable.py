from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from hailtop.batch.job import Job

from lrs_annotation.scripts.svs import vcf_to_annotate_cohort_mt


def annotate_cohort_jobs_svs(
    vcf_path: Path,
    out_mt_path: Path,
    gencode_gtf_path: Path,
    checkpoint_prefix: Path,
    job_attrs: dict | None = None,
) -> Job:
    """
    Annotate cohort VCF for seqr loader, SNPs and Indels.
    """
    j = get_batch().new_job('Annotate cohort', job_attrs)
    j.image(config_retrieve(['workflow', 'driver_image']))
    gencode_gz = get_batch().read_input(gencode_gtf_path)
    j.command(
        f"""
        python3 {vcf_to_annotate_cohort_mt.__file__} \\
            --vcf_path {vcf_path} \\
            --out_mt_path {out_mt_path} \\
            --gencode_gz {gencode_gz} \\
            --checkpoint_prefix {checkpoint_prefix} \
        """
    )
    return j
