"""
Jobs to screen a LongTR VCF against STRchive disease loci
and produce standalone HTML reports, optionally subset by loci list.
"""

from hailtop.batch.job import Job

from cpg_utils import Path, config, hail_batch

from lrs_annotation.scripts import longtr_pathogenic


def longtr_pathogenic_report(
    vcf_path: str,
    outputs: dict[str, Path],
    job_attrs: dict[str, str],
    loci_lists: dict[str, list[str]] | None = None,
) -> Job:
    """
    Run the LongTR pathogenic screening script on a VCF file.

    If loci_lists is provided, produces one HTML+JSON per report type.
    Otherwise produces a single report with all loci.
    """
    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_job('LongTR Pathogenic Report', job_attrs | {'tool': 'longtr_pathogenic'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.storage('10GB')

    local_vcf = batch_instance.read_input(vcf_path)
    strchive_json = batch_instance.read_input(config.config_retrieve(['references', 'strchive_json']))
    longtr_bed = batch_instance.read_input(config.config_retrieve(['references', 'strchive_longtr_bed']))

    if not loci_lists:
        job.command(f"""
    python3 {longtr_pathogenic.__file__} \\
        --vcf_path {local_vcf} \\
        --strchive_json {strchive_json} \\
        --longtr_bed {longtr_bed} \\
        --output_html {job.html} \\
        --output_json {job.json}
    """)
        batch_instance.write_output(job.html, str(outputs['html']))
        batch_instance.write_output(job.json, str(outputs['json']))
    else:
        for list_name, loci in loci_lists.items():
            html_rg = job[f'{list_name}_html']
            json_rg = job[f'{list_name}_json']
            loci_str = ' '.join(loci)
            job.command(f"""
    python3 {longtr_pathogenic.__file__} \\
        --vcf_path {local_vcf} \\
        --strchive_json {strchive_json} \\
        --longtr_bed {longtr_bed} \\
        --output_html {html_rg} \\
        --output_json {json_rg} \\
        --report_type {list_name} \\
        --loci_list {loci_str}
    """)
            batch_instance.write_output(html_rg, str(outputs[f'{list_name}_html']))
            batch_instance.write_output(json_rg, str(outputs[f'{list_name}_json']))

    return job
