from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import Batch

from lrs_annotation.utils import get_resource_overrides_for_job

def reformat_snps_indels_vcf_with_bcftools(
    batch: Batch,
    job_name: str,
    job_attrs: dict,
    vcf_path: str,
    lrs_sg_id_mapping_path: str,
):
    """
    Reformat SNPs and Indels VCFs using bcftools.
    """
    # Input files
    local_vcf = batch.read_input(vcf_path)
    local_id_mapping = batch.read_input(lrs_sg_id_mapping_path)

    # Required file for normalisation
    ref_fasta = config_retrieve(['workflow', 'ref_fasta'])
    fasta = batch.read_input_group(**{'fa': ref_fasta, 'fa.fai': f'{ref_fasta}.fai'})['fa']

    # Use BCFtools to reheader the VCF, replacing the LRS IDs with the SG IDs
    reformatting_job = batch.new_job(job_name, job_attrs)
    reformatting_job.declare_resource_group(
        vcf_out={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}
    )
    reformatting_job.image(image=config_retrieve['images', 'bcftools'])

    resource_overrides = get_resource_overrides_for_job('reformat_snps_indels_vcf')

    reformatting_job.storage(resource_overrides.get('storage_gb', '10Gi'))
    reformatting_job.command(
        f'bcftools view -Ov {local_vcf} | bcftools reheader --samples {local_id_mapping} '
        f'-o {reformatting_job.reheadered}',
    )
    # Normalise, sort, bgzip, and index the VCF
    reformatting_job.command(
        f'bcftools norm -m -any -f {fasta} -c s {reformatting_job.reheadered} | '
        f'bcftools sort | bgzip -c > {reformatting_job.vcf_out["vcf.bgz"]}',
    )
    reformatting_job.command(f'tabix {reformatting_job.vcf_out["vcf.bgz"]}')

    return reformatting_job
