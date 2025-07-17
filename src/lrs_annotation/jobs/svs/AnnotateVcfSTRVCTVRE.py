
import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import Path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch
from lrs_annotation.utils import get_references


def annotate_strvctvre_job(
    input_vcf: hb.ResourceFile,
    output_path: Path,
    job_attrs: dict | None = None,
    name: str = 'AnnotateVcfWithStrvctvre',
) -> Job:
    """

    Args:
        input_vcf (ResourceFile): part of a resource group with the corresponding index
        output_path ():
        job_attrs (dict|None): job attributes
        name (str): name of the job

    Returns:
        The Strvctvre job
    """

    job_attrs = job_attrs or {}
    strv_job = get_batch().new_job('StrVCTVRE', job_attrs | {'tool': 'strvctvre'})

    strv_job.image(image_path('strvctvre'))
    strv_job.storage(config_retrieve(['resource_overrides', name, 'storage'], '10Gi'))
    strv_job.memory(config_retrieve(['resource_overrides', name, 'memory'], '16Gi'))

    strvctvre_phylop = get_references(['strvctvre_phylop'])['strvctvre_phylop']
    assert isinstance(strvctvre_phylop, str)

    local_phylop = get_batch().read_input(strvctvre_phylop)

    strv_job.declare_resource_group(output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})

    # run strvctvre
    strv_job.command(
        f'python StrVCTVRE.py -i {input_vcf} -o {strv_job.output["vcf.gz"]} -f vcf -p {local_phylop}',
    )
    strv_job.command(f'tabix {strv_job.output["vcf.gz"]}')

    get_batch().write_output(strv_job.output, str(output_path).replace('.vcf.gz', ''))
    return strv_job
