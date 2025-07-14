
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from lrs_annotation.scripts.svs import sniffles_vcf_modifier

def sniffles_modify_vcf(
    vcf_path: Path,
    ref_fa_path: Path,
    sex_mapping_file_path: Path,
    job_name: str,
    job_attrs: dict | None = None,
) -> Job:
    """
    Call the sniffles VCF modifier script to modify the SVs VCF file before reformatting.
    """
    j = get_batch().new_job(job_name, job_attrs)

    vcf = get_batch().read_input(vcf_path)
    fasta = get_batch().read_input_group(**{'fa': ref_fa_path, 'fa.fai': f'{ref_fa_path}.fai'})['fa']
    sex_mapping = get_batch().read_input(sex_mapping_file_path)
    j.declare_resource_group(
        vcf_out={'vcf.gz': '{root}.vcf.gz'}
    )

    j.image(config_retrieve(['workflow', 'driver_image']))
    j.storage('10Gi')
    j.command(
        f"""
        python3 {sniffles_vcf_modifier.__file__} \\
            --vcf_path {vcf} \\
            --vcf_out {j.vcf_out["vcf.gz"]} \\
            --fa {fasta} \\
            --sex_mapping_file {sex_mapping}
        """
    )
    return j
