"""
Job to screen a LongTR VCF against STRchive disease loci
and produce a standalone HTML report.
"""

from hailtop.batch.job import Job

from cpg_utils import Path, config, hail_batch

from lrs_annotation.scripts import longtr_pathogenic


def longtr_pathogenic_report(
    vcf_path: str,
    output: Path,
    job_attrs: dict[str, str],
) -> Job:
    """
    Run the LongTR pathogenic screening script on a VCF file.

    Localizes the VCF and STRchive reference files onto the VM,
    then runs the script to produce an HTML report.
    """
    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_job('LongTR Pathogenic Report', job_attrs | {'tool': 'longtr_pathogenic'})

    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.storage('10GB')

    local_vcf = batch_instance.read_input(vcf_path)
    strchive_json = batch_instance.read_input(config.config_retrieve(['references', 'strchive_json']))
    longtr_bed = batch_instance.read_input(config.config_retrieve(['references', 'strchive_longtr_bed']))

    job.command(f"""
    python3 {longtr_pathogenic.__file__} \\
        --vcf_path {local_vcf} \\
        --strchive_json {strchive_json} \\
        --longtr_bed {longtr_bed} \\
        --output {job.output_file}
    """)

    batch_instance.write_output(job.output_file, str(output))
    return job
