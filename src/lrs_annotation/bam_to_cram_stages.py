"""
Workflow for annotating long-read SVs data into a seqr-ready format.
"""

from cpg_flow import stage, targets
from cpg_flow.filetypes import CramPath
from cpg_utils import Path
from cpg_utils.config import config_retrieve, reference_path
from cpg_utils.hail_batch import get_batch

from lrs_annotation.jobs import BamToCram, SomalierExtract





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
        cram_path = sg.cram or sg.dataset.prefix() / 'cram' / f'{sg.id}.cram'
        return {'cram': cram_path.path, 'crai': cram_path.index_path}

    def queue_jobs(self, sg: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Using the existing `bam_to_cram` function from the `jobs` module.
        """
        input_bam = get_batch().read_input_group(bam=str(sg.alignment_input))
        job = BamToCram.bam_to_cram(
            b=get_batch(),
            input_bam=input_bam,
            job_attrs=self.get_job_attrs(sg),
        )
        get_batch().write_output(job.sorted_cram, str(self.expected_outputs(sg)['cram']).removesuffix('.cram'))

        return self.make_outputs(sg, data=self.expected_outputs(sg), jobs=[job])


@stage.stage(required_stages=[ConvertBamToCram])
class CramQcSomalier(stage.SequencingGroupStage):
    """Run somalier extract on a CRAM file."""

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> Path:
        return sg.cram or sg.dataset.prefix() / 'cram' / f'{sg.id}.cram.somalier'

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(sequencing_group)

        cram = inputs.as_str(sequencing_group, ConvertBamToCram, 'cram')

        jobs = SomalierExtract.extract_somalier(
            cram_path=cram,
            output=output,
            job_attrs=self.get_job_attrs(sequencing_group),
        )

        return self.make_outputs(sequencing_group, data=output, jobs=jobs)
