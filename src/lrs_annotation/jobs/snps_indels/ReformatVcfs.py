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

    Reheader the VCF to replace the LRS IDs with the SG IDs with bcftools reheader.
    Normalise the VCF with the reference genome and bcftools norm.
    Finally, sort and index the VCF with bcftools sort.
    """
    # Input files
    local_vcf = batch.read_input(vcf_path)
    local_id_mapping = batch.read_input(lrs_sg_id_mapping_path)

    # Required file for normalisation
    ref_fasta = config_retrieve(['workflow', 'ref_fasta'])
    fasta = batch.read_input_group(**{'fa': ref_fasta, 'fa.fai': f'{ref_fasta}.fai'})['fa']

    # Use BCFtools to reheader the VCF, replacing the LRS IDs with the SG IDs
    reformatting_job = batch.new_job(job_name, job_attrs)
    reformatting_job.image(image=config_retrieve(['images', 'bcftools']))
    reformatting_job.declare_resource_group(
        vcf_out={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    reformatting_job = get_resource_overrides_for_job(reformatting_job, 'reformat_snps_indels_vcf')
    reformatting_job.command(
        f'bcftools reheader --samples {local_id_mapping} {local_vcf} | '
        f'bcftools norm --check-refs s -m -any -f {fasta} -c s -Ou | '
        f'bcftools sort -Oz -W=tbi - -o {reformatting_job.vcf_out["vcf.gz"]}'
    )
    return reformatting_job
