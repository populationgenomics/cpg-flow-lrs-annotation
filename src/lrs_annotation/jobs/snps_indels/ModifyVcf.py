from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from hailtop.batch.job import Job

from lrs_annotation.utils import get_resource_overrides_for_job


def bcftools_reformat(
    vcf_path: str,
    job_name: str,
    job_attrs: dict,
    lrs_sg_id_mapping_path: str,
) -> Job:
    """
    Reformat SNPs and Indels VCFs using bcftools.

    First validates that sample IDs in the VCF exist in the mapping file.
    Reheader the VCF to replace the LRS IDs with the SG IDs with bcftools reheader.
    Normalise the VCF with the reference genome and bcftools norm.
    Uppercases the allele strings with awk.
    Finally, sort and index the VCF with bcftools sort.
    """
    # Input files
    local_vcf = get_batch().read_input(vcf_path)
    local_id_mapping = get_batch().read_input(lrs_sg_id_mapping_path)

    # Required file for normalisation
    ref_fasta = config_retrieve(['workflow', 'ref_fasta'])
    fasta = get_batch().read_input_group(**{'fa': ref_fasta, 'fa.fai': f'{ref_fasta}.fai'})['fa']

    # Use BCFtools to reheader the VCF, replacing the LRS IDs with the SG IDs
    job = get_batch().new_job(job_name, job_attrs)
    job.image(image=config_retrieve(['images', 'bcftools']))
    job.declare_resource_group(vcf_out={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
    job = get_resource_overrides_for_job(job, 'reformat_snps_indels_vcf')

    # First, validate that sample ID(s) from the input VCF exist in the reheadering / mapping file
    job.command(f"""
        # Extract sample IDs from VCF
        SAMPLE_IDS=$(bcftools query -l {local_vcf})

        # Check each sample ID exists in mapping file
        for sample_id in $SAMPLE_IDS; do
            if ! grep -q "^$sample_id\\s" {local_id_mapping}; then
                echo "ERROR: Sample ID '$sample_id' not found in mapping file {lrs_sg_id_mapping_path}"
                exit 1
            else
                echo "✓ Sample ID '$sample_id' found in mapping file"
            fi
        done

        echo "All sample IDs validated successfully"
    """)

    job.command(
        f'bcftools reheader --samples {local_id_mapping} {local_vcf} | '
        f'bcftools norm -m -any -f {fasta} -c s -Ou | '
        f'bcftools sort -Ov - | '
        # Uppercase the allele strings with awk
        'awk \'BEGIN{OFS="\t"} /^#/ {print; next} {$4=toupper($4); $5=toupper($5); print}\' | '
        f'bgzip > {job.vcf_out["vcf.gz"]} && '
        f'tabix -p vcf {job.vcf_out["vcf.gz"]}'
    )
    return job
