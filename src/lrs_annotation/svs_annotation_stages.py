"""
Workflow for annotating long-read SVs data into a seqr-ready format.
"""

from cpg_flow import stage, targets, workflow
from cpg_flow.utils import tshirt_mt_sizing
from cpg_flow.workflow import get_multicohort
from cpg_utils import Path, to_path
from cpg_utils.config import AR_GUID_NAME, config_retrieve, try_get_ar_guid
from cpg_utils.hail_batch import get_batch
from google.api_core.exceptions import PermissionDenied
from inputs import get_sgs_from_datasets, query_for_lrs_mappings, query_for_lrs_vcfs
from jobs.ExportMtToElasticsearch import export_mt_to_elasticsearch
from jobs.svs.AnnotateCohortMatrixtable import annotate_cohort_jobs_svs
from jobs.svs.AnnotateDatasetMatrixtable import annotate_dataset_jobs_sv
from jobs.svs.AnnotateVcfGatk import queue_annotate_sv_jobs
from jobs.svs.AnnotateVcfSTRVCTVRE import annotate_strvctvre_job
from jobs.svs.MergeVcfs import merge_svs_vcf_with_bcftools
from jobs.svs.ModifySvVcf import modify_sv_vcf
from jobs.svs.ReformatVcfs import reformat_svs_vcf_with_bcftools
from jobs.svs.WriteCleanedPedFile import make_clean_combined_ped
from loguru import logger
from utils import (
    es_password,
    get_dataset_name,
    get_dataset_names,
    get_query_filter_from_config,
    write_mapping_to_file,
)


@stage.stage
class WriteLrsIdToSgAndSexMappingFiles(stage.MultiCohortStage):
    """
    Write the LRS ID to SG ID and LRS ID to sex mappings to files.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'lrs_sg_id_mapping': self.prefix / 'lrs_sg_id_mapping.txt',
            'lrs_id_sex_mapping': self.prefix / 'lrs_id_sex_mapping.txt',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Write the LRS ID to SG ID mapping to a file, and the LRS ID to sex mapping to another file.
        The first is used by bcftools reheader to update the sample IDs in the VCFs, the second is used by the
        SV VCF modifier script when calculating copy number.
        """
        outputs = self.expected_outputs(multicohort)
        sgid_mapping_file_path = outputs['lrs_sg_id_mapping']
        sex_mapping_file_path = outputs['lrs_id_sex_mapping']

        lrs_mapping = query_for_lrs_mappings(
            dataset_names=get_dataset_names([d.name for d in multicohort.get_datasets()]),
            sequencing_types=get_query_filter_from_config('sequencing_types', make_tuple=False),
            sequencing_platforms=get_query_filter_from_config('sequencing_platforms', make_tuple=False),
        )
        lrs_sg_id_mapping = {lrs_id: mapping['sg_id'] for lrs_id, mapping in lrs_mapping.items()}
        lrs_sex_mapping = {lrs_id: mapping['sex'] for lrs_id, mapping in lrs_mapping.items()}

        logger.info(f'Writing LRS ID to SG ID mapping to {sgid_mapping_file_path}')
        write_mapping_to_file(lrs_sg_id_mapping, sgid_mapping_file_path)
        logger.info(f'Writing LRS ID to participant sex mapping to {sex_mapping_file_path}')
        write_mapping_to_file(lrs_sex_mapping, sex_mapping_file_path)

        return self.make_outputs(multicohort, data=self.expected_outputs(multicohort))


@stage.stage
class WriteCleanedPedFile(stage.MultiCohortStage):
    """
    Write a cleaned PED file for the cohort, concatenating all samples across all datasets with the reference panel.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'ped_file': self.prefix / 'ped_with_ref_panel.ped',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Create a cleaned PED file for the cohort.
        """
        outputs = self.expected_outputs(multicohort)
        logger.info(f'Writing cleaned PED file to {outputs["ped_file"]}')
        make_clean_combined_ped(multicohort.write_ped_file(), outputs['ped_file'])
        return self.make_outputs(multicohort, data=outputs)


@stage.stage(required_stages=[WriteLrsIdToSgAndSexMappingFiles])
class ModifySVsVcf(stage.SequencingGroupStage):
    """
    Modify the long-read SV VCFs, saving them to temp storage for use in the next stage.
    - This is explicitly skipped for the parents in trio joint-called VCFs
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
        sgid_prefix = sequencing_group.dataset.tmp_prefix() / 'svs' / 'modified_vcfs'
        return {
            'vcf': sgid_prefix / f'{sequencing_group.id}_modified.vcf.gz',
        }

    def queue_jobs(self, sg: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Uses a python job to change the VCF contents, making the following changes:
        - Replace REF strings with the reference base
        - Replace ALT strings with symbolic "<TYPE>" derived from the SVTYPE INFO field
        - Add a CN field to the FORMAT field, filled with copy number values based on the genotype and sex
        - Adds a unique ID to each record for compatibility with GATK SV sorting
        """
        multicohort_datasets = [ds.name for ds in get_multicohort().get_datasets()]
        sg_ids = []
        sg_vcfs = {}
        for ds in multicohort_datasets:
            sgs = query_for_lrs_vcfs(dataset_name=ds)
            sg_ids.extend(sgs['sg_ids'])
            sg_vcfs.update(sgs['vcfs'])

        assert not set(get_multicohort().get_sequencing_group_ids()) - set(sg_ids), (
            'SGs in the multicohort do not have VCFs matching the filter criteria: '
            f'{set(get_multicohort().get_sequencing_group_ids()) - set(sg_ids)}'
            ' Adjust the query filter or the input cohorts'
        )

        if sg.id not in sg_vcfs:
            return None

        vcf_path: str = sg_vcfs[sg.id]['vcf']
        lrs_id_sex_map = inputs.as_path(get_multicohort(), WriteLrsIdToSgAndSexMappingFiles, 'lrs_id_sex_mapping')

        joint_called = sg_vcfs[sg.id]['meta'].get('joint_called', False)
        job_name = f'Modify {"joint-called " if joint_called else ""}{vcf_path} prior to reformatting'

        mod_job = modify_sv_vcf(
            vcf_path=to_path(vcf_path),
            ref_fa_path=config_retrieve(['workflow', 'ref_fasta']),
            sex_mapping_file_path=lrs_id_sex_map,
            job_name=job_name,
            job_attrs={'tool': 'python'},
        )

        # write from temp storage into GCP
        expected_outputs = self.expected_outputs(sg)
        get_batch().write_output(mod_job.vcf_out, str(expected_outputs['vcf']).removesuffix('.vcf.gz'))

        return self.make_outputs(target=sg, jobs=[mod_job], data=expected_outputs)


@stage.stage(required_stages=[ModifySVsVcf, WriteLrsIdToSgAndSexMappingFiles])
class ReformatSVsVcfWithBcftools(stage.SequencingGroupStage):
    """
    Reformat the modified long-read SV VCFs using BCFtools
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
        sgid_prefix = sequencing_group.dataset.tmp_prefix() / 'svs' / 'reformatted_vcfs'
        return {
            'vcf': sgid_prefix / f'{sequencing_group.id}_reformatted.vcf.gz',
            'index': sgid_prefix / f'{sequencing_group.id}_reformatted.vcf.gz.tbi',
        }

    def queue_jobs(self, sg: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        A python job to change the VCF contents
        - Use bcftools job to reheader the VCF with the replacement sample IDs, normalise it, and then sort
        - Then block-gzip and index it
        """
        sg_vcfs = query_for_lrs_vcfs(dataset_name=get_dataset_name(sg.dataset.name))['vcfs']
        if sg.id not in sg_vcfs:
            return None

        # Input VCF and reheadering file
        vcf_path = inputs.as_path(sg, ModifySVsVcf, 'vcf')
        lrs_sg_id_mapping = inputs.as_path(get_multicohort(), WriteLrsIdToSgAndSexMappingFiles, 'lrs_sg_id_mapping')

        reformatting_job = reformat_svs_vcf_with_bcftools(
            batch=get_batch(),
            job_name=f'Reformat SVs VCF for {sg.id}',
            job_attrs={'tool': 'bcftools'},
            vcf_path=vcf_path,
            lrs_sg_id_mapping_path=lrs_sg_id_mapping,
        )

        outputs = self.expected_outputs(sg)
        get_batch().write_output(reformatting_job.vcf_out, str(outputs['vcf']).removesuffix('.vcf.gz'))

        return self.make_outputs(target=sg, jobs=[reformatting_job], data=outputs)


@stage.stage(required_stages=ReformatSVsVcfWithBcftools)
class MergeSVsVcfsWithBcftools(stage.MultiCohortStage):
    """
    Merge the reformatted SVs VCFs together with bcftools
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'vcf': self.tmp_prefix / 'svs' / 'merged_reformatted.vcf.gz',
            'index': self.tmp_prefix / 'svs' / 'merged_reformatted.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Use bcftools to merge all the VCFs, and then fill in the tags (requires bcftools 1.18+)
        """
        sgs = get_sgs_from_datasets([d.name for d in multicohort.get_datasets()])
        # Get the reformatted VCFs from the previous stage
        reformatted_vcfs = inputs.as_dict_by_target(ReformatSVsVcfWithBcftools)
        reformatted_vcfs = {sg_id: vcf for sg_id, vcf in reformatted_vcfs.items() if sg_id in sgs['vcfs']}

        if len(reformatted_vcfs) == 1:
            logger.info('Only one VCF found, skipping merge')
            return None

        vcf_paths = [
            str(reformatted_vcfs.get(sgid, {}).get('vcf'))
            for sgid in multicohort.get_sequencing_group_ids()
            if sgid in reformatted_vcfs
        ]

        merge_job = merge_svs_vcf_with_bcftools(batch=get_batch(), vcf_paths=vcf_paths, job_attrs={'tool': 'bcftools'})

        outputs = self.expected_outputs(multicohort)
        get_batch().write_output(merge_job.output, str(outputs['vcf']).removesuffix('.vcf.gz'))

        return self.make_outputs(multicohort, data=outputs, jobs=merge_job)


@stage.stage(required_stages=[WriteCleanedPedFile, ModifySVsVcf, ReformatSVsVcfWithBcftools, MergeSVsVcfsWithBcftools])
class AnnotateSVsWithGatk(stage.MultiCohortStage):
    """
    Annotates the merged long-read SVs VCF with GATK-SV.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'annotated_vcf': self.tmp_prefix / 'svs' / 'gatk_annotated.vcf.gz',
            'annotated_vcf_index': self.tmp_prefix / 'svs' / 'gatk_annotated.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Use the GATK-SV wrapper to annotate this merged VCF.
        """
        outputs = self.expected_outputs(multicohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        if len(inputs.as_dict_by_target(ReformatSVsVcfWithBcftools)) == 1:
            sg = multicohort.get_sequencing_groups()[0]
            input_vcf = inputs.as_path(sg, ReformatSVsVcfWithBcftools, 'vcf')
        else:
            input_vcf = inputs.as_path(multicohort, MergeSVsVcfsWithBcftools, 'vcf')

        jobs = queue_annotate_sv_jobs(
            dataset_name=multicohort.analysis_dataset,
            multicohort_name=multicohort.name,
            multicohort_ped_file_path=inputs.as_path(multicohort, WriteCleanedPedFile, 'ped_file'),
            input_vcf=input_vcf,
            outputs=outputs,
            labels=billing_labels,
        )
        return self.make_outputs(multicohort, jobs=jobs, data=outputs)


@stage.stage(required_stages=AnnotateSVsWithGatk)
class AnnotateSVsWithStrvctvre(stage.MultiCohortStage):
    """
    Annotate the long-read SVs VCF with STRVCTVRE.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'strvctvre_vcf': self.tmp_prefix / 'svs' / 'strvctvre_annotated.vcf.gz',
            'strvctvre_vcf_index': self.tmp_prefix / 'svs' / 'strvctvre_annotated.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Queue a job to annotate the VCF from the previous stage with STRVCTVRE.
        """
        outputs = self.expected_outputs(multicohort)
        input_dict = inputs.as_dict(multicohort, AnnotateSVsWithGatk)
        input_vcf = get_batch().read_input_group(
            vcf=str(input_dict['annotated_vcf']),
            vcf_index=str(input_dict['annotated_vcf_index']),
        )['vcf']

        strvctvre_job = annotate_strvctvre_job(
            input_vcf=input_vcf,
            output_path=outputs['strvctvre_vcf'],
            job_attrs=self.get_job_attrs(),
        )

        return self.make_outputs(multicohort, jobs=strvctvre_job, data=outputs)


@stage.stage(required_stages=AnnotateSVsWithStrvctvre)
class AnnotateCohortSVsMtFromVcfWithHail(stage.MultiCohortStage):
    """
    First step to transform annotated SV callset data into a seqr ready format
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        """
        Expected to write a matrix table.
        """
        return {'mt': self.tmp_prefix / 'svs' / 'cohort.mt'}

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        queue job(s) to rearrange the annotations prior to Seqr transformation
        """
        outputs = self.expected_outputs(multicohort)
        vcf_path = inputs.as_path(target=multicohort, stage=AnnotateSVsWithStrvctvre, key='strvctvre_vcf')

        job = annotate_cohort_jobs_svs(
            vcf_path=vcf_path,
            out_mt_path=outputs['mt'],
            gencode_gtf_path=config_retrieve(['workflow', 'gencode_gtf_file']),
            checkpoint_prefix=self.tmp_prefix / 'svs' / 'checkpoints',
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage.stage(required_stages=[AnnotateCohortSVsMtFromVcfWithHail], analysis_type='sv', analysis_keys=['mt'])
class SubsetSVsMtToDatasetWithHail(stage.DatasetStage):
    """
    Subset the MT to be this Dataset only
    Then work up all the genotype values
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict:
        """
        Expected to generate a matrix table
        """
        sg_hash = workflow.get_workflow().output_version
        return {'mt': (dataset.prefix() / 'mt' / f'LongReadSV-{sg_hash}-{dataset.name}.mt')}

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Subsets the whole MT to this cohort only
        Then brings a range of genotype data into row annotations

        Args:
            dataset (Dataset): SGIDs specific to this dataset/project
            inputs ():
        """
        # only create matrixtables for datasets specified in the config
        # or, if this config option is not set, run for all datasets
        eligible_datasets = config_retrieve(['workflow', 'write_mt_for_datasets'], default=None)
        if eligible_datasets is not None and dataset.name not in eligible_datasets:
            logger.info(f'Skipping MT writing for {dataset}')
            return self.make_outputs(dataset)
        mt_path = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortSVsMtFromVcfWithHail, key='mt')
        outputs = self.expected_outputs(dataset)

        sg_hash = workflow.get_workflow().output_version
        checkpoint_prefix = dataset.tmp_prefix() / sg_hash / 'svs' / 'mt' / 'checkpoints'

        jobs = annotate_dataset_jobs_sv(
            mt_path=mt_path,
            sg_ids=dataset.get_sequencing_group_ids(),
            out_mt_path=outputs['mt'],
            tmp_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage.stage(
    required_stages=[SubsetSVsMtToDatasetWithHail],
    analysis_type='es-index',
    analysis_keys=['done_flag'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'SV'},  # noqa: ARG005
)
class ExportSVsMtToElasticIndex(stage.DatasetStage):
    """
    Create a Seqr index
    https://github.com/populationgenomics/metamist/issues/539
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sg_hash = workflow.get_workflow().output_version
        sequencing_type = config_retrieve(['workflow', 'sequencing_type'])
        index_name = f'{dataset.name}-{sequencing_type}-LR-SV-{sg_hash}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'svs' / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Uses the non-DataProc MT-to-ES conversion script
        """
        # only create the elasticsearch index for the datasets specified in the config
        # or, if this config option is not set, create it for all datasets
        eligible_datasets = config_retrieve(['workflow', 'create_es_index_for_datasets'], default=None)
        if eligible_datasets is not None and dataset.name not in eligible_datasets:
            logger.info(f'Skipping ES index creation for {dataset}')
            return None
        # try to generate a password here - we'll find out inside the script anyway, but
        # by that point we'd already have localised the MT, wasting time and money
        try:
            es_password()
        except PermissionDenied:
            logger.warning(f'No permission to access ES password, skipping for {dataset}')
            return self.make_outputs(dataset)
        except KeyError:
            logger.warning(f'ES section not in config, skipping for {dataset}')
            return self.make_outputs(dataset)

        outputs = self.expected_outputs(dataset)

        # get the absolute path to the MT
        mt_path = str(inputs.as_path(target=dataset, stage=SubsetSVsMtToDatasetWithHail, key='mt'))

        # get the expected outputs as Strings
        index_name = str(outputs['index_name'])
        flag_name = str(outputs['done_flag'])
        # and just the name, used after localisation
        mt_name = mt_path.split('/')[-1]

        req_storage = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(get_sgs_from_datasets([dataset.name])),
        )

        job = export_mt_to_elasticsearch(
            batch=get_batch(),
            mt_path=mt_path,
            index_name=index_name,
            flag_name=flag_name,
            req_storage=req_storage,
            job_name=f'Export {index_name} from {mt_name}',
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs, jobs=job)
