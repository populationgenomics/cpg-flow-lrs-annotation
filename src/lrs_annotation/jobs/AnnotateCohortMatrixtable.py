
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from lrs_annotation.scripts import annotate_cohort_snps_indels

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
    init_batch_args: dict[str, str | int] = {}
    annotate_cohort_workflow = config_retrieve(['workflow', 'annotate_cohort'], {})

    # Memory parameters
    for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
        if annotate_cohort_workflow.get(config_key):
            init_batch_args[batch_key] = 'highmem'
    # Cores parameter
    for key in ['driver_cores', 'worker_cores']:
        if annotate_cohort_workflow.get(key):
            init_batch_args[key] = annotate_cohort_workflow[key]

    j = get_batch().new_job('Annotate cohort', job_attrs)
    j.image(config_retrieve(['workflow', 'driver_image']))
    # TODO: find a way to get the driver / worker resources into this from the config, like query_command did
    j.command(
        f"""
        python3 {annotate_cohort_snps_indels.__file__} \\
            --vcf_path {vcf_path} \\
            --out_mt_path {out_mt_path} \\
            --vep_ht_path {vep_ht_path} \\
            --checkpoint_prefix {checkpoint_prefix} \\
            --remove_invalid_contigs
        """
    )
    return j
