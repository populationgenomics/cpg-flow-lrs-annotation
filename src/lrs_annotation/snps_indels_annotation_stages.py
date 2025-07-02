"""
This file exists to define all the Stages for the workflow.
The logic for each stage can be contained here (if it is not too complex),
or can be delegated to a separate file in jobs.

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

Workflow for finding long-read SNPs_Indels files, updating file contents, merging, and annotating the results
"""

from loguru import logger

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch

from cpg_flow import stage, targets
from cpg_flow.utils import tshirt_mt_sizing
from cpg_flow.workflow import (
    get_multicohort,
    get_workflow,
)

from jobs.SplitMergedVcfAndGetSitesOnlyForVep import split_merged_vcf_and_get_sitesonly_vcfs_for_vep
from jobs.AnnotateWithVep import add_vep_jobs
from jobs.AnnotateCohortMatrixtable import annotate_cohort_jobs_snps_indels
from jobs.AnnotateDatasetMatrixtable import annotate_dataset_jobs

from utils import (
    get_config_options_as_tuple,
    query_for_lrs_vcfs,
    query_for_lrs_to_sg_id_and_sex_mapping,
    write_mapping_to_file,
    joint_calling_scatter_count,
    es_password,
)

from google.api_core.exceptions import PermissionDenied

@stage.stage
class WriteLrsIdToSgIdMappingFile(stage.MultiCohortStage):
    """
    Write the LRS ID to SG ID mapping to a file
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'lrs_sgid_mapping': self.prefix / 'lrs_sgid_mapping.txt',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Write the LRS ID to SG ID mapping to a file
        This is used by bcftools reheader to update the sample IDs in the VCFs
        """
        output = self.expected_outputs(multicohort)

        lrs_sgid_mapping, _ = query_for_lrs_to_sg_id_and_sex_mapping([d.name for d in multicohort.get_datasets()])
        mapping_file_path = self.prefix / 'lrs_sgid_mapping.txt'
        logger.info(f'Writing LRS ID to SG ID mapping to {mapping_file_path}')
        write_mapping_to_file(lrs_sgid_mapping, mapping_file_path)

        return self.make_outputs(multicohort, data=output)


@stage.stage(required_stages=[WriteLrsIdToSgIdMappingFile])
class ReformatSnpsIndelsVcfWithBcftools(stage.SequencingGroupStage):
    """
    Take each of the long-read SNPs Indels VCFs, and re-format the contents
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
        sgid_prefix = sequencing_group.dataset.prefix() / 'pacbio' / 'modified_vcfs'
        return {
            'vcf': sgid_prefix / f'{sequencing_group.id}_reformatted_snps_indels.vcf.bgz',
            'index': sgid_prefix / f'{sequencing_group.id}_reformatted_snps_indels.vcf.bgz.tbi',
        }

    def queue_jobs(self, sg: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        A python job to change the VCF contents
        - Use bcftools job to reheader the VCF with the replacement sample IDs, normalise it, and then sort
        - Then block-gzip and index it
        - This is explicitly skipped for the parents in trio joint-called VCFs
        """
        sg_vcfs = query_for_lrs_vcfs(
            dataset_name=sg.dataset.name,
            sequencing_types=get_config_options_as_tuple('sequencing_types'),
            variant_types=get_config_options_as_tuple('variant_types'),
            variant_callers=get_config_options_as_tuple('variant_callers'),
            pipeface_versions=get_config_options_as_tuple('pipeface_versions'),
            joint_called=config_retrieve(['workflow', 'lrs_annotation', 'joint_called'], default=False),
            verbose=config_retrieve(['workflow', 'lrs_annotation', 'verbose'], default=False),
        )
        if sg.id not in sg_vcfs:
            return None

        joint_called = sg_vcfs[sg.id]['meta'].get('joint_called', False)

        # Input VCF
        sg_vcf: str = sg_vcfs[sg.id]['output']
        local_vcf = get_batch().read_input(sg_vcf)

        # Required file for reheadering
        lrs_sg_id_mapping = inputs.as_path(get_multicohort(), WriteLrsIdToSgIdMappingFile, 'lrs_sgid_mapping')
        local_id_mapping = get_batch().read_input(lrs_sg_id_mapping)

        # Required file for normalisation
        ref_fasta = config_retrieve(['workflow', 'ref_fasta'])
        fasta = get_batch().read_input_group(**{'fa': ref_fasta, 'fa.fai': f'{ref_fasta}.fai'})['fa']

        # Use BCFtools to reheader the VCF, replacing the LRS IDs with the SG IDs
        reformatting_job = get_batch().new_job(
            f'Reheadering, Normalising, BGZipping, and Indexing for {sg.id}: '
            f'{"joint-called " if joint_called else ""}{sg_vcf}',
            {'tool': 'bcftools'},
        )
        reformatting_job.declare_resource_group(
            vcf_out={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}
        )
        reformatting_job.image(image=image_path('bcftools'))
        reformatting_job.storage('10Gi')
        reformatting_job.command(
            f'bcftools view -Ov {local_vcf} | bcftools reheader --samples {local_id_mapping} '
            f'-o {reformatting_job.reheadered}',
        )
        # Normalise, sort, bgzip, and index the VCF
        reformatting_job.command(
            f'bcftools norm -m -any -f {fasta} -c s {reformatting_job.reheadered} | '
            f'bcftools sort | bgzip -c > {reformatting_job.vcf_out["vcf.bgz"]}',
        )
        reformatting_job.command(f'tabix {reformatting_job.vcf_out["vcf.bgz"]}')

        outputs = self.expected_outputs(sg)
        get_batch().write_output(reformatting_job.vcf_out, str(outputs['vcf']).removesuffix('.vcf.bgz'))

        return self.make_outputs(target=sg, jobs=[reformatting_job], data=outputs)


@stage.stage(required_stages=ReformatSnpsIndelsVcfWithBcftools)
class MergeReformattedSnpsIndelsVcfsWithBcftools(stage.MultiCohortStage):
    """
    Find all the reformatted VCFs and merge them together with bcftools
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'vcf': self.prefix / 'merged_reformatted_snps_indels.vcf.bgz',
            'index': self.prefix / 'merged_reformatted_snps_indels.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Use bcftools to merge all the VCFs, and then fill in the tags (requires bcftools 1.18+)
        """
        # Get the reformatted VCFs
        reformatted_vcfs = inputs.as_dict_by_target(ReformatSnpsIndelsVcfWithBcftools)
        if len(reformatted_vcfs) == 1:
            logger.info('Only one VCF found, skipping merge')
            return None

        batch_vcfs = [
            get_batch().read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz']
            for each_vcf in [
                str(reformatted_vcfs.get(sgid, {}).get('vcf'))
                for sgid in multicohort.get_sequencing_group_ids()
                if sgid in reformatted_vcfs
            ]
        ]

        merge_job = get_batch().new_job('Merge Long-Read SNPs Indels calls', attributes={'tool': 'bcftools'})
        merge_job.image(image=image_path('bcftools_121'))

        merge_job_resources = config_retrieve(
            ['workflow', 'lrs_annotation', 'resource_overrides', 'merge_job'],
            default={
                'cpu': 4,
                'memory': '16Gi',
                'storage': '50Gi',
            }
        )
        merge_job.cpu(merge_job_resources['cpu'])
        merge_job.memory(merge_job_resources['memory'])
        merge_job.storage(merge_job_resources['storage'])
        merge_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # BCFtools options breakdown:
        #   --threads: number of threads to use
        #   -m: merge strategy (none means multiple records for multiallelic sites)
        #   -0: assume genotypes at missing sites are 0/0
        #   -Oz: bgzip output (compressed VCF)
        #   -o: output file name
        #   --write-index: write index file (only for bcftools 1.18+. Occasionally bugged for < 1.21)
        #   +fill-tags: plugin to compute and fill in the INFO tags (AF, AN, AC)
        merge_job.command(
            f'bcftools merge {" ".join(batch_vcfs)} --threads 4 -m none -0 | '
            f'bcftools +fill-tags -Oz -o {merge_job.output["vcf.bgz"]} --write-index=tbi -- -t AF,AN,AC'
        )

        outputs = self.expected_outputs(multicohort)
        get_batch().write_output(merge_job.output, str(outputs['vcf']).removesuffix('.vcf.bgz'))

        return self.make_outputs(multicohort, data=outputs, jobs=merge_job)


@stage.stage(required_stages=[ReformatSnpsIndelsVcfWithBcftools, MergeReformattedSnpsIndelsVcfsWithBcftools])
class SplitMergedVcfIntoSitesOnlyVcfsWithGatk(stage.MultiCohortStage):
    """
    Get the site-only VCFs from the merged VCF by splitting the VCF with GATK SelectVariants,
    then extracting the sites-only VCFs for VEP annotation.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict:
        """
        Generate site-only VCFs from the merged VCF.
        """
        return {
            'siteonly': to_path(self.prefix / 'siteonly.vcf.gz'),
            'siteonly_part_pattern': str(self.prefix / 'siteonly_parts' / 'part{idx}.vcf.gz'),
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Submit jobs.
        """
        jobs = []
        outputs = self.expected_outputs(multicohort)
        scatter_count = joint_calling_scatter_count(len(multicohort.get_sequencing_groups()))
        out_siteonly_vcf_part_paths = [
            to_path(outputs['siteonly_part_pattern'].format(idx=idx)) for idx in range(scatter_count)
        ]

        intervals_path = None
        if config_retrieve(['workflow', 'intervals_path'], default=None):
            intervals_path = to_path(config_retrieve(['workflow', 'intervals_path']))

        exclude_intervals_path = None
        if config_retrieve(['workflow', 'exclude_intervals_path'], default=None):
            exclude_intervals_path = to_path(config_retrieve(['workflow', 'exclude_intervals_path']))

        if len(inputs.as_dict_by_target(ReformatSnpsIndelsVcfWithBcftools)) == 1:
            merged_vcf_path = inputs.as_path(multicohort, ReformatSnpsIndelsVcfWithBcftools, 'vcf')
        else:
            merged_vcf_path = inputs.as_path(multicohort, MergeReformattedSnpsIndelsVcfsWithBcftools, 'vcf')

        vcf_jobs = split_merged_vcf_and_get_sitesonly_vcfs_for_vep(
            b=get_batch(),
            scatter_count=scatter_count,
            merged_vcf_path=merged_vcf_path,
            tmp_bucket=self.tmp_prefix / 'tmp',
            out_siteonly_vcf_part_paths=out_siteonly_vcf_part_paths,
            intervals_path=intervals_path,
            exclude_intervals_path=exclude_intervals_path,
            job_attrs=self.get_job_attrs(),
        )
        jobs.extend(vcf_jobs)

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[SplitMergedVcfIntoSitesOnlyVcfsWithGatk])
class VepLongReadAnnotation(stage.MultiCohortStage):
    """
    Run VEP on the long-read site-only VCFs and write out a Hail table.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort):
        """
        Expected to write a hail table.
        """
        return {'ht': self.prefix / 'long_read' / 'vep.ht'}

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Submit jobs.
        """
        outputs = self.expected_outputs(multicohort)
        scatter_count = joint_calling_scatter_count(len(multicohort.get_sequencing_groups()))
        input_siteonly_vcf_part_paths = [
            to_path(
                inputs.as_str(
                    stage=SplitMergedVcfIntoSitesOnlyVcfsWithGatk,
                    target=multicohort,
                    key='siteonly_part_pattern',
                ).format(idx=idx),
            )
            for idx in range(scatter_count)
        ]

        jobs = add_vep_jobs(
            get_batch(),
            input_vcfs=input_siteonly_vcf_part_paths,
            out_path=outputs['ht'],
            tmp_prefix=self.tmp_prefix / 'tmp',
            job_attrs=self.get_job_attrs(),
            scatter_count=scatter_count,
        )
        return self.make_outputs(multicohort, outputs, jobs)


@stage.stage(required_stages=[
    ReformatSnpsIndelsVcfWithBcftools,
    MergeReformattedSnpsIndelsVcfsWithBcftools,
    VepLongReadAnnotation]
)
class AnnotateCohortLongReadSnpsIndelsWithHail(stage.MultiCohortStage):
    """
    First step to transform annotated SNPs Indels callset data into a seqr ready format
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        """
        Expected to write a matrix table.
        """
        return {'mt': self.prefix / 'cohort_snps_indels.mt'}

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        queue job(s) to rearrange the annotations prior to Seqr transformation
        """
        outputs = self.expected_outputs(multicohort)

        if len(inputs.as_dict_by_target(ReformatSnpsIndelsVcfWithBcftools)) == 1:
            vcf_path = inputs.as_path(multicohort, ReformatSnpsIndelsVcfWithBcftools, 'vcf')
        else:
            vcf_path = inputs.as_path(target=multicohort, stage=MergeReformattedSnpsIndelsVcfsWithBcftools, key='vcf')

        vep_ht_path = inputs.as_path(target=multicohort, stage=VepLongReadAnnotation, key='ht')

        job = annotate_cohort_jobs_snps_indels(
            vcf_path=vcf_path,
            out_mt_path=outputs['mt'],
            vep_ht_path=vep_ht_path,
            checkpoint_prefix=self.tmp_prefix / 'checkpoints',
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage.stage(required_stages=[AnnotateCohortLongReadSnpsIndelsWithHail], analysis_type='custom', analysis_keys=['mt'])
class SubsetAnnotateCohortLongReadSnpsIndelsMtToDatasetWithHail(stage.DatasetStage):
    """
    Subset the MT to be this Dataset only
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict:
        """
        Expected to generate a matrix table
        """

        return {
            'mt': (dataset.prefix() / 'mt' / f'LongReadSNPsIndels-{get_workflow().output_version}-{dataset.name}.mt'),
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Subsets the whole MT to this cohort only
        Then brings a range of genotype data into row annotations

        Args:
            dataset (Dataset): SGIDs specific to this dataset/project
            inputs ():
        """
        mt_path = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortLongReadSnpsIndelsWithHail, key='mt')

        outputs = self.expected_outputs(dataset)

        checkpoint_prefix = dataset.tmp_prefix() / dataset.name / 'checkpoints'

        jobs = annotate_dataset_jobs(
            mt_path=mt_path,
            sequencing_group_ids=dataset.get_sequencing_group_ids(),
            out_mt_path=outputs['mt'],
            tmp_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage.stage(
    required_stages=[SubsetAnnotateCohortLongReadSnpsIndelsMtToDatasetWithHail],
    analysis_type='es-index',
    analysis_keys=['index_name'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'SNV_INDEL'},  # noqa: ARG005
)
class ExportLongReadSnpsIndelsMtToElasticIndex(stage.DatasetStage):
    """
    Create a Seqr index
    https://github.com/populationgenomics/metamist/issues/539
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = config_retrieve(['workflow', 'sequencing_type'])
        index_name = f'{dataset.name}-{sequencing_type}-LR-SNPsIndels-{get_workflow().run_timestamp}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Uses the non-DataProc MT-to-ES conversion script
        """

        # try to generate a password here - we'll find out inside the script anyway, but
        # by that point we'd already have localised the MT, wasting time and money
        try:
            _es_password_string = es_password()
        except PermissionDenied:
            logger.warning(f'No permission to access ES password, skipping for {dataset}')
            return self.make_outputs(dataset)
        except KeyError:
            logger.warning(f'ES section not in config, skipping for {dataset}')
            return self.make_outputs(dataset)

        outputs = self.expected_outputs(dataset)

        # get the absolute path to the MT
        mt_path = str(
            inputs.as_path(target=dataset, stage=SubsetAnnotateCohortLongReadSnpsIndelsMtToDatasetWithHail, key='mt')
        )
        # and just the name, used after localisation
        mt_name = mt_path.split('/')[-1]

        # get the expected outputs as Strings
        index_name = str(outputs['index_name'])
        flag_name = str(outputs['done_flag'])

        job = get_batch().new_job(f'Generate {index_name} from {mt_path}')
        req_storage = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(dataset.get_sequencing_group_ids()),
        )
        # set all job attributes in one bash
        job.cpu(4).memory('lowmem').storage(f'{req_storage}Gi').image(config_retrieve(['workflow', 'driver_image']))

        # localise the MT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')

        # run the export from the localised MT - this job writes no new data, just transforms and exports over network
        job.command(f'mt_to_es --mt_path "${{BATCH_TMPDIR}}/{mt_name}" --index {index_name} --flag {flag_name}')

        return self.make_outputs(dataset, data=outputs['index_name'], jobs=job)
