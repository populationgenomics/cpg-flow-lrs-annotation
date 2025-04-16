"""
trivially simple python job example, using a utility constant
"""

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from workflow_name.python_jobs.PrintPreviousJobOutputInAPythonJob import print_file_contents

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def set_up_printing_python_job(input_file: str, output_file: str) -> 'Job':
    """
    This is a simple example of a job that prints the contents of a file.
    This is the logic for the stage, calling the pythonJob to do the work

    Args:
        input_file (str): the path to the file to print
        output_file (str): the path to write the result to
    """

    # localise the file
    local_input = get_batch().read_input(input_file)

    # run the PythonJob
    job = get_batch().new_python_job(f'Read {input_file}')
    job.image(config_retrieve(['workflow', 'driver_image']))
    pyjob_output = job.call(
        print_file_contents,
        local_input,
    )
    get_batch().write_output(pyjob_output.as_str(), output_file)
    return job
