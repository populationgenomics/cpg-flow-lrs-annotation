
from hailtop.batch.job import Job
from cpg_utils import Path
from cpg_utils.config import config_retrieve
from lrs_annotation.utils import add_gatk_sv_jobs, get_references, get_images
from typing import Any

def queue_annotate_sv_jobs(
    dataset_name: str,
    multicohort_name: str,
    multicohort_ped_file_path: Path,
    input_vcf: Path,
    outputs: dict,
    labels: dict[str, str] | None = None,
) -> list[Job] | Job | None:
    """
    Helper function to queue jobs for SV annotation
    Enables common access to the same Annotation WDL for CNV & SV
    """
    input_dict: dict[str, Any] = {
        'vcf': input_vcf,
        'prefix': multicohort_name,
        'ped_file': multicohort_ped_file_path,
        'sv_per_shard': 5000,
        'external_af_population': config_retrieve(['references', 'gatk_sv', 'external_af_population']),
        'external_af_ref_prefix': config_retrieve(['references', 'gatk_sv', 'external_af_ref_bed_prefix']),
        'external_af_ref_bed': config_retrieve(['references', 'gnomad_sv']),
    }

    input_dict |= get_references(
        [
            'noncoding_bed',
            'protein_coding_gtf',
            {'contig_list': 'primary_contigs_list'},
        ],
    )

    # images!
    input_dict |= get_images(['sv_pipeline_docker', 'sv_base_mini_docker', 'gatk_docker'])
    return add_gatk_sv_jobs(
        dataset_name=dataset_name,
        wfl_name='AnnotateVcf',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels=labels,
    )
