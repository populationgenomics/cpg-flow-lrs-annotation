"""
Convert BAM to CRAM.
"""

from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, image_path, reference_path
from cpg_utils.hail_batch import Batch, command
from cpg_flow.resources import STANDARD


def bam_to_cram(
    b: Batch,
    input_bam: ResourceGroup,
    job_attrs: dict | None = None,
) -> tuple[Job, ResourceGroup]:
    """
    Convert a BAM file to a CRAM file.
    """
    assert isinstance(input_bam, ResourceGroup)

    job_name = 'bam_to_cram'

    convert_tool = 'samtools_view'
    j_attrs = (job_attrs or {}) | {'label': job_name, 'tool': convert_tool}
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    reference_fasta_path = to_path(config_retrieve(['workflow', 'ref_fasta'], reference_path('broad/ref_fasta')))

    # Get fasta file
    fasta = b.read_input_group(
        fasta=reference_fasta_path,
        fasta_fai=f'{reference_fasta_path}.fai',
    )

    # Set resource requirements
    res = STANDARD.set_resources(
        j=j,
        ncpu=config_retrieve(['resource_overrides', 'bam_to_cram', 'cpu_cores'], 1),
        storage_gb=config_retrieve(['workflow', 'resource_overrides', 'bam_to_cram', 'storage_gib'], 50),
    )

    j.declare_resource_group(
        sorted_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        },
    )

    cmd = f"""
    samtools view -@ {res.get_nthreads() - 1} -T {fasta.fasta} -C {input_bam.bam} | \
    tee {j.sorted_cram["cram"]} | samtools index -@ {res.get_nthreads() - 1} - {j.sorted_cram["cram.crai"]}
    """
    j.command(command(cmd, monitor_space=True))

    return j, j.sorted_cram
