
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from lrs_annotation.scripts import annotate_cohort_snps_indels
from lrs_annotation.utils import get_init_batch_args_for_job

def annotate_cohort_jobs_snps_indels(
    vcf_path: Path,
    out_mt_path: Path,
    vep_ht_path: Path,
    checkpoint_prefix: Path,
    job_attrs: dict | None = None,
) -> Job:
    """
    Annotate cohort VCF for seqr loader, SNPs and Indels.
    """
    init_batch_args = get_init_batch_args_for_job('annotate_cohort_snps_indels')

    j = get_batch().new_job('Annotate cohort', job_attrs)
    j.image(config_retrieve(['workflow', 'driver_image']))
    j.command(
        f"""
        python3 {annotate_cohort_snps_indels.__file__} \\
            --vcf_path {vcf_path} \\
            --out_mt_path {out_mt_path} \\
            --vep_ht_path {vep_ht_path} \\
            --checkpoint_prefix {checkpoint_prefix} \\
            --remove_invalid_contigs \\
        """
        f'  --init_batch_args {init_batch_args}' if init_batch_args else ' \\'

    )
    return j
