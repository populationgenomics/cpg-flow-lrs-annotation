"""
Job for somalier self-relatedness check.
Runs somalier extract on each input file (using cached .somalier where available),
then somalier relate to compare all fingerprints.
"""

import re

from cpg_flow.utils import exists
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import Job

_CRAM_PATTERN = re.compile(r'\.cram$')


def somalier_self_check(
    sg_files: dict[str, str],
    somalier_dir: Path,
    output_prefix: Path,
    job_attrs: dict[str, str],
    sg_id: str,
    participant_id: str,
    dataset_name: str,
    ped_file: str | None = None,
) -> list[Job]:
    """
    Run somalier extract + relate as a single job.

    For each SG, uses cpg_flow.utils.exists() at orchestration time to check
    if a .somalier file already exists in GCS. If so, it's loaded via read_input
    instead of re-extracting. For SGs that need extraction, runs somalier extract.
    Assumes single-sample input files (VCF, gVCF, or CRAM).
    Then runs somalier relate on all .somalier files.

    Args:
        sg_files: mapping of SG ID -> input file path (VCF, gVCF, or CRAM)
        somalier_dir: GCS directory for .somalier output files
        output_prefix: GCS path prefix for relate output files
        job_attrs: job attributes for Hail Batch
        ped_file: optional PED file for pedigree-aware relatedness (future use)
    """
    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_job('Somalier self-check', job_attrs | {'tool': 'somalier'})
    job.image(config.config_retrieve(['images', 'somalier']))

    has_cram = any(_CRAM_PATTERN.search(path) for path in sg_files.values())
    if has_cram:
        storage_gb = config.config_retrieve(
            ['workflow', 'resource_overrides', 'bam_to_cram', 'storage_gib'],
            50,
        )
        job.storage(f'{storage_gb}GB')
    else:
        job.storage('10GB')

    ref = hail_batch.fasta_res_group(batch_instance)
    sites = batch_instance.read_input(config.config_retrieve(['references', 'somalier_sites']))

    # For each SG, check at orchestration time if .somalier already exists.
    # If so, use read_input to load it. Otherwise, localise the input file for extraction.
    cached_somalier: dict[str, str] = {}  # sg_id -> localised cached .somalier path
    needs_extract: dict[str, str] = {}  # sg_id -> localised input file path

    for sg_id, file_path in sg_files.items():
        somalier_output_path = somalier_dir / f'{sg_id}.somalier'
        if exists(somalier_output_path):
            cached_somalier[sg_id] = str(batch_instance.read_input(str(somalier_output_path)))
        else:
            is_cram = bool(_CRAM_PATTERN.search(file_path))
            if is_cram:
                localised = batch_instance.read_input_group(
                    cram=file_path,
                    crai=f'{file_path}.crai',
                ).cram
            else:
                localised = batch_instance.read_input_group(
                    **{'vcf.gz': file_path, 'vcf.gz.tbi': f'{file_path}.tbi'},
                )['vcf.gz']
            needs_extract[sg_id] = str(localised)

    # Build shell commands: copy cached files, extract new ones
    cmds = ['mkdir -p extracted/']

    for sg_id, cached_path in cached_somalier.items():
        cmds.append(f'cp {cached_path} extracted/{sg_id}.somalier')

    for sg_id, localised_path in needs_extract.items():
        cmds.append(
            f'export SOMALIER_SAMPLE_NAME={sg_id}\n'
            f'mkdir -p tmp_{sg_id}\n'
            f'somalier extract -d tmp_{sg_id}/ --sites {sites} -f {ref.base} {localised_path}\n'
            f'mv tmp_{sg_id}/*.somalier extracted/{sg_id}.somalier\n'
            f'cp extracted/{sg_id}.somalier {job[f"somalier_{sg_id}"]}'
        )

    ped_arg = ''
    if ped_file:
        localised_ped = batch_instance.read_input(ped_file)
        ped_arg = f'--ped {localised_ped}'

    cmds.extend(
        [
            '',
            '# Relate all extracted fingerprints',
            f'somalier relate {ped_arg} -o results extracted/*.somalier',
            f'mv results.pairs.tsv {job.pairs_tsv}',
            f'mv results.samples.tsv {job.samples_tsv}',
            f'mv results.html {job.html}',
        ]
    )

    job.command('\n'.join(cmds))

    for sg_id in needs_extract:
        batch_instance.write_output(job[f'somalier_{sg_id}'], str(somalier_dir / f'{sg_id}.somalier'))

    # Write relate outputs
    batch_instance.write_output(job.pairs_tsv, str(output_prefix) + '.pairs.tsv')
    batch_instance.write_output(job.samples_tsv, str(output_prefix) + '.samples.tsv')
    batch_instance.write_output(job.html, str(output_prefix) + '.html')

    # Check self-relatedness and alert via Slack if kinship < threshold
    kinship_threshold = config.config_retrieve(
        ['workflow', 'somalier_self_check', 'kinship_threshold'], 0.9,
    )
    check_job = batch_instance.new_job('Somalier self-check alert', job_attrs)
    check_job.image(config.config_retrieve(['workflow', 'driver_image']))
    check_job.depends_on(job)

    hail_batch.copy_common_env(check_job)
    hail_batch.authenticate_cloud_credentials_in_job(check_job)
    check_job.command(f"""\
    python3 -m lrs_annotation.scripts.check_self_relatedness \\
    --pairs-tsv {job.pairs_tsv} \\
    --sg-id {sg_id} \\
    --participant-id {participant_id} \\
    --dataset {dataset_name} \\
    --kinship-threshold {kinship_threshold}
    """)

    return [job, check_job]
