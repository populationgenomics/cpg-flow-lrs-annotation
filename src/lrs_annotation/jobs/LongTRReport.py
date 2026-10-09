"""
Jobs to screen a LongTR VCF against STRchive disease loci
and produce standalone HTML reports, optionally subset by loci list.
"""

from hailtop.batch.job import Job

from cpg_utils import Path, config, hail_batch


def longtr_pathogenic_report(
    vcf_path: str,
    outputs: dict[str, Path],
    job_attrs: dict[str, str],
    loci_lists: dict[str, list[str]],
    sample_id: str,
    *,
    sex: str = 'unknown',
    birth_year: str = '',
    age_of_onset: str = '',
    hpo_terms: str = '',
) -> Job:
    """
    Run the LongTR pathogenic screening script on a VCF file, one HTML+JSON per loci list.
    """
    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_job('LongTR Pathogenic Report', job_attrs)
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.storage('10GB')

    # the script fetches each locus by region, so the tabix index has to be localised next to
    # the VCF - the group keys carry the extensions so htslib finds it at <vcf>.tbi
    local_vcf = batch_instance.read_input_group(
        **{'vcf.gz': vcf_path, 'vcf.gz.tbi': f'{vcf_path}.tbi'},
    )['vcf.gz']
    strchive_json = batch_instance.read_input(config.config_retrieve(['references', 'strchive_json']))
    custom_loci = batch_instance.read_input(config.config_retrieve(['references', 'strchive_custom_json']))

    # phenotype fields are optional per participant, so only pass what metamist holds
    phenotype_args = []
    if birth_year:
        phenotype_args.extend(['--birth_year', f"'{birth_year}'"])
    if age_of_onset:
        phenotype_args.extend(['--age_of_onset', f"'{age_of_onset}'"])
    if hpo_terms:
        phenotype_args.extend(['--hpo_terms', f"'{hpo_terms}'"])
    phenotype_arg = ' '.join(phenotype_args)

    for list_name, loci in loci_lists.items():
        html_rg = job[f'{list_name}_html']
        json_rg = job[f'{list_name}_json']
        loci_str = ' '.join(loci)
        # --loci_list is nargs='+', so every other argument has to precede it
        job.command(f"""
    python3 -m lrs_annotation.scripts.longtr_pathogenic \\
        --vcf_path {local_vcf} \\
        --strchive_json {strchive_json} \\
        --custom_loci_json {custom_loci} \\
        --output_html {html_rg} \\
        --output_json {json_rg} \\
        --report_type {list_name} \\
        --sex {sex} \\
        --sample_id '{sample_id}' {phenotype_arg} \\
        --loci_list {loci_str}
    """)
        batch_instance.write_output(html_rg, str(outputs[f'{list_name}_html']))
        batch_instance.write_output(json_rg, str(outputs[f'{list_name}_json']))

    return job
