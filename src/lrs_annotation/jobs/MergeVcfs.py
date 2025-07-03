from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import Batch

from lrs_annotation.utils import get_resource_overrides_for_job

def merge_snps_indels_vcf_with_bcftools(
    batch: Batch,
    vcf_paths: str,
):
    """
    Reformat SNPs and Indels VCFs using bcftools.
    """
    # Input files
    batch_vcfs = [
        batch.read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz']
        for each_vcf in vcf_paths
    ]

    merge_job = batch.new_job('Merge Long-Read SNPs Indels calls', attributes={'tool': 'bcftools'})
    merge_job.image(image=config_retrieve(['images', 'bcftools_121']))

    resource_overrides = get_resource_overrides_for_job('merge_snps_indels_vcfs')

    merge_job.cpu(resource_overrides.get('cpu', 4))
    merge_job.memory(resource_overrides.get('memory', '16Gi'))
    merge_job.storage(resource_overrides.get('storage', '50Gi'))
    merge_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    # BCFtools options breakdown:
    #   --threads: number of threads to use
    #   -m: merge strategy (none means multiple records for multiallelic sites)
    #   -0: assume genotypes at missing sites are 0/0
    #   -Oz: bgzip output (compressed VCF)
    #   -o: output file name
    #   --write-index: write index file (only for bcftools 1.18+. Occasionally bugged for < 1.21)
    #   +fill-tags: plugin to compute and fill in the INFO tags (AF, AN, AC)
    merge_job.command(
        f'bcftools merge {" ".join(batch_vcfs)} --threads 4 -m none -0 | '
        f'bcftools +fill-tags -Oz -o {merge_job.output["vcf.bgz"]} --write-index=tbi -- -t AF,AN,AC'
    )

    return merge_job
