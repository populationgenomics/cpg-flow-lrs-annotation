"""
Workflow for annotating long-read SVs data into a seqr-ready format.
"""

from cpg_flow import stage, targets
from cpg_flow.filetypes import CramPath
from cpg_utils import Path
from cpg_utils.config import config_retrieve, reference_path
from cpg_utils.hail_batch import get_batch
from jobs.BamToCram import bam_to_cram


def make_long_read_cram_path(sg: targets.SequencingGroup) -> CramPath:
    """
    Path to a CRAM file. Not checking its existence here.
    """
    path: Path = sg.dataset.prefix() / 'long_read' / 'cram' / f'{sg.id}.cram'
    return CramPath(
        path=path,
        index_path=path.with_suffix('.cram.crai'),
        reference_assembly=config_retrieve(['workflow', 'ref_fasta'], reference_path('broad/ref_fasta')),
    )


@stage.stage(
    analysis_type=config_retrieve(['workflow', 'bam_to_cram', 'analysis_type'], 'cram'),
    analysis_keys=['cram'],
)
class ConvertBamToCram(stage.SequencingGroupStage):
    """
    Convert a PacBio BAM to a CRAM file.
    """

    def expected_outputs(self, sg: targets.SequencingGroup) -> dict[str, Path]:
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        cram_path = sg.cram or make_long_read_cram_path(sg)
        return {'cram': cram_path.path, 'crai': cram_path.index_path}

    def queue_jobs(self, sg: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Using the existing `bam_to_cram` function from the `jobs` module.
        """
        input_bam = get_batch().read_input_group(bam=str(sg.alignment_input))
        job = bam_to_cram(
            b=get_batch(),
            input_bam=input_bam,
            job_attrs=self.get_job_attrs(sg),
        )
        get_batch().write_output(job.sorted_cram, str(self.expected_outputs(sg)['cram']).removesuffix('.cram'))

        return self.make_outputs(sg, data=self.expected_outputs(sg), jobs=[job])
