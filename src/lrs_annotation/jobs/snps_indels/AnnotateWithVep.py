
import hailtop.batch as hb
from hailtop.batch.job import Job

from textwrap import dedent
from typing import Literal
from lrs_annotation.scripts.snps_indels import vep_json_to_ht
from lrs_annotation.jobs.snps_indels.GatherVcfs import gather_vcfs
from cpg_flow.utils import can_reuse
from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, reference_path


def add_vep_jobs(
    b: hb.Batch,
    input_vcfs: list[Path],
    tmp_prefix: Path,
    scatter_count: int,
    out_path: Path | None = None,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Runs VEP on provided VCF. Writes a VCF into `out_path` by default,
    unless `out_path` ends with ".ht", in which case writes a Hail table.
    """
    to_hail_table = out_path and out_path.suffix == '.ht'
    if not to_hail_table:
        assert str(out_path).endswith('.vcf.gz'), out_path

    if out_path and can_reuse(out_path):
        return []

    jobs: list[Job] = []

    if len(input_vcfs) == 0:
        raise ValueError('No input VCFs provided')

    # read all input VCFs as resource groups
    input_vcf_resources: list[hb.ResourceGroup] = [
        b.read_input_group(**{'vcf.gz': str(path), 'vcf.gz.tbi': str(path) + '.tbi'}) for path in input_vcfs
    ]

    result_parts_bucket = tmp_prefix / 'vep' / 'parts'
    result_part_paths = []
    for idx, resource in enumerate(input_vcf_resources):
        if to_hail_table:
            result_part_path = result_parts_bucket / f'part{idx + 1}.jsonl'
        else:
            result_part_path = result_parts_bucket / f'part{idx + 1}.vcf.gz'
        result_part_paths.append(result_part_path)
        if can_reuse(result_part_path):
            continue

        # noinspection PyTypeChecker
        vep_one_job = vep_one(
            b,
            vcf=resource['vcf.gz'],
            out_format='json' if to_hail_table else 'vcf',
            out_path=result_part_paths[idx],
            job_attrs=(job_attrs or {}) | {'part': f'{idx + 1}/{scatter_count}'},
        )
        if vep_one_job:
            jobs.append(vep_one_job)

    if to_hail_table:
        j = gather_vep_json_to_ht(
            b=b,
            vep_results_paths=result_part_paths,
            out_path=out_path,
            job_attrs=job_attrs,
        )
        gather_jobs = [j]
    elif scatter_count != 1:
        assert len(result_part_paths) == scatter_count
        gather_jobs = gather_vcfs(b=b, input_vcfs=result_part_paths, out_vcf_path=out_path)
    else:
        print('no need to merge VEP results')
        gather_jobs = []

    for j in gather_jobs:
        j.depends_on(*jobs)

    jobs.extend(gather_jobs)

    return jobs


def gather_vep_json_to_ht(
    b: hb.Batch,
    vep_results_paths: list[Path],
    out_path: Path,
    job_attrs: dict | None = None,
) -> Job:
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table using a Batch job.
    """
    j = b.new_job('VEP', job_attrs)
    j.image(config_retrieve(['workflow', 'driver_image']))
    j.command(
        f"""
        python3 {vep_json_to_ht.__file__} \\
        --vep_results_paths {' '.join(map(str, vep_results_paths))} \
        --out_path {out_path}
        """
    )
    return j


def vep_one(
    b: hb.Batch,
    vcf: Path | hb.ResourceFile,
    out_path: Path | None = None,
    out_format: Literal['vcf', 'json'] = 'vcf',
    job_attrs: dict | None = None,
) -> Job | None:
    """
    Run a single VEP job.
    """
    if out_path and can_reuse(out_path):
        return None

    # check that the cache and image for VEP 110 exist
    vep_image = config_retrieve(['images', 'vep_110'])
    vep_mount_path = to_path(reference_path('vep_110_mount'))
    assert all([vep_image, vep_mount_path])

    j = b.new_job('VEP', (job_attrs or {}) | {'tool': 'VEP 110'})
    j.image(vep_image)

    # vep is single threaded, with a middling memory requirement
    # during test it can exceed 8GB, so we'll give it 16GB
    j.memory('16Gi').storage('15Gi').cpu(1)

    if not isinstance(vcf, hb.ResourceFile):
        vcf = b.read_input(str(vcf))

    if out_format == 'vcf':
        j.declare_resource_group(output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
        assert isinstance(j.output, hb.ResourceGroup)
        output = j.output['vcf.gz']
    else:
        assert isinstance(j.output, hb.ResourceFile)
        output = j.output

    # gcsfuse works only with the root bucket, without prefix:
    data_mount = to_path(f'/{vep_mount_path.drive}')
    j.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
    vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])

    # assume for now that only VEP 105 has a non-standard install location
    loftee_conf = {
        'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
        'conservation_file': f'{vep_dir}/loftee.sql',
        'loftee_path': '$VEP_DIR_PLUGINS',
    }

    # plugin only present in VEP 110 build
    alpha_missense_plugin = f'--plugin AlphaMissense,file={vep_dir}/AlphaMissense_hg38.tsv.gz '

    # VCF annotation doesn't utilise the aggregated Seqr reference data, including spliceAI
    # SpliceAI requires both indel and SNV files to be present (~100GB), untested
    use_splice_ai = config_retrieve(['workflow', 'spliceai_plugin'], False)
    vcf_plugins = (
        (
            f'--plugin SpliceAI,snv={vep_dir}/spliceai_scores.raw.snv.hg38.vcf.gz,'
            f'indel={vep_dir}/spliceai_scores.raw.indel.hg38.vcf.gz '
        )
        if (use_splice_ai and out_format == 'vcf')
        else ''
    )

    cmd = f"""
    set -x
    vep \\
    --format vcf \\
    --{out_format} {'--compress_output bgzip' if out_format == 'vcf' else ''} \\
    -o {output} \\
    -i {vcf} \\
    --everything \\
    --mane_select \\
    --allele_number \\
    --minimal \\
    --species homo_sapiens \\
    --cache --offline --assembly GRCh38 \\
    --dir_cache {vep_dir}/vep/ \\
    --fasta /cpg-common-main/references/vep/110/mount/vep/homo_sapiens/110/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \\
    {alpha_missense_plugin} \\
    --plugin LoF,{','.join(f'{k}:{v}' for k, v in loftee_conf.items())} \\
    --plugin UTRAnnotator,file=$UTR38 {vcf_plugins}
    """

    if out_format == 'vcf':
        cmd += f'tabix -p vcf {output}'

    j.command(dedent(cmd))

    if out_path:
        b.write_output(j.output, str(out_path).replace('.vcf.gz', ''))

    return j
