"""
This file exists to define all the Stages for the workflow.
The logic for each stage should be delegated to a separate file in jobs.

In this repository template I've segregated stages (this file), jobs (logic) and python_jobs (methods to be invoked only
as python jobs). This means that even the most trivial Stage has it's logic in a separate file. In reality Stages are
not likely to consist of single line print statements, so it's good to enforces strict separation of logic.

- A Stage contains the definition of inputs, outputs, and it's place in the workflow DAG (via required_stages)
- A job encapsulates all the logic required to complete the Stage, and is contained in a separate file
- A stage will call its corresponding job to queu operations, and the job will return the Hail Batch Job object(s)

Naming conventions for Stages are not enforced, but a series of recommendations have been made here:

https://cpg-populationanalysis.atlassian.net/wiki/spaces/ST/pages/185597962/Pipeline+Naming+Convention+Specification

A suggested naming convention for a stages is:
  - PascalCase (each word capitalized, no hyphens or underscores)
  - If the phrase contains an initialism (e.g. VCF), only the first character should be capitalised
  - Verb + Subject (noun) + Preposition + Direct Object (noun)  TODO(anyone): please correct my grammar is this is false
  e.g. AlignShortReadsWithBowtie2, or MakeSitesOnlyVcfWithBcftools
  - This becomes self-explanatory when reading the code and output folders

Each Stage should be a Class, and should inherit from one of
  - SequencingGroupStage
  - DatasetStage
  - CohortStage
  - MultiCohortStage
"""

from typing import TYPE_CHECKING

from workflow_name.jobs.DoSomethingGenericWithBash import echo_statement_to_file
from workflow_name.jobs.PrintPreviousJobOutputInAPythonJob import set_up_printing_python_job

from cpg_flow.stage import MultiCohortStage, stage

if TYPE_CHECKING:
    # Path is a classic return type for a Stage, and is a shortcut for [CloudPath | pathlib.Path]
    from cpg_utils import Path
    from cpg_flow.targets import MultiCohort, StageInput, StageOutput


@stage()
class DoSomethingGenericWithBash(MultiCohortStage):
    """
    This is a generic stage that runs a bash command.
    """

    def expected_outputs(self, multicohort: 'MultiCohort') -> 'Path':
        """
        This is where we define the expected outputs for this stage.
        """
        # self.prefix() is a more concise shortcut for multicohort.analysis_dataset_bucket/ StageName / Hash
        return multicohort.analysis_dataset.prefix(category='tmp') / self.name / 'output.txt'

    def queue_jobs(self, multicohort: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':  # noqa: ARG002
        """
        This is where we generate jobs for this stage.
        """

        # locate the intended output path
        outputs = self.expected_outputs(multicohort)

        # generate the output
        j = echo_statement_to_file('Hello World!', str(outputs))

        # return the jobs and outputs
        return self.make_outputs(multicohort, data=outputs, jobs=j)


@stage(required_stages=[DoSomethingGenericWithBash])
class PrintPreviousJobOutputInAPythonJob(MultiCohortStage):
    """
    This is a stage that cats the output of a previous stage to the logs.
    It is implemented here as a Stage which calls a job method - that job then creates and queues a PythonJob
    """

    def expected_outputs(self, multicohort: 'MultiCohort') -> 'Path':
        return multicohort.analysis_dataset.prefix(category='tmp') / self.name / 'cat.txt'

    def queue_jobs(self, multicohort: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        # get the previous stage's output
        previous_output = inputs.as_str(multicohort, DoSomethingGenericWithBash)

        # generate the expected output path
        outputs = self.expected_outputs(multicohort)

        # run the PythonJob
        job = set_up_printing_python_job(
            input_file=previous_output,
            output_file=str(outputs),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=job)
