"""
Workflow for annotating long-read SNPs and Indels data into a seqr-ready format.
"""

from loguru import logger

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from cpg_flow import stage, targets, workflow
from cpg_flow.utils import tshirt_mt_sizing
from cpg_flow.workflow import get_multicohort

from jobs.snps_indels.AnnotateCohortMatrixtable import annotate_cohort_jobs_snps_indels
from jobs.snps_indels.AnnotateDatasetMatrixtable import annotate_dataset_jobs
from jobs.snps_indels.AnnotateWithVep import add_vep_jobs
from jobs.snps_indels.MergeVcfs import merge_snps_indels_vcf_with_bcftools
from jobs.snps_indels.ModifyVcf import bcftools_reformat
from jobs.snps_indels.SplitMergedVcfAndGetSitesOnlyForVep import split_merged_vcf_and_get_sitesonly_vcfs_for_vep
from jobs.ExportMtToElasticsearch import export_mt_to_elasticsearch
from inputs import (
    get_sgs_from_datasets,
    query_for_lrs_vcfs,
    query_for_lrs_mappings
)

from utils import (
    get_dataset_names,
    get_family_sequencing_groups,
    get_query_filter_from_config,
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
            'lrs_sg_id_mapping': self.prefix / 'lrs_sg_id_mapping.txt',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Write the LRS ID to SG ID mapping to a file
        This is used by bcftools reheader to update the sample IDs in the VCFs
        """
        outputs = self.expected_outputs(multicohort)
        mapping_file_path = outputs['lrs_sg_id_mapping']

        lrs_mapping = query_for_lrs_mappings(
            dataset_names=get_dataset_names([d.name for d in multicohort.get_datasets()]),
            sequencing_types=get_query_filter_from_config('sequencing_types', make_tuple=False),
            sequencing_platforms=get_query_filter_from_config('sequencing_platforms', make_tuple=False)
        )
        lrs_sg_id_mapping = {lrs_id: mapping['sg_id'] for lrs_id, mapping in lrs_mapping.items()}

        logger.info(f'Writing LRS ID to SG ID mapping to {mapping_file_path}')
        write_mapping_to_file(lrs_sg_id_mapping, mapping_file_path)

        return self.make_outputs(multicohort, data=self.expected_outputs(multicohort))


@stage.stage(required_stages=[WriteLrsIdToSgIdMappingFile])
class ModifyVcf(stage.SequencingGroupStage):
    """
    Modify the long-read SNPs Indels VCFs as a pre-processing step before merging.
    """
    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
        sgid_prefix = sequencing_group.dataset.tmp_prefix() / 'snps_indels' / 'reformatted_vcfs'
        return {
            'vcf': sgid_prefix / f'{sequencing_group.id}_reformatted.vcf.gz',
            'index': sgid_prefix / f'{sequencing_group.id}_reformatted.vcf.gz.tbi',
        }

    def queue_jobs(self, sg: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        """
        - Use bcftools job to reheader the VCF with the replacement sample IDs, normalise it, and then sort
        - Then block-gzip and index it
        - This is explicitly skipped for the parents in trio joint-called VCFs
        """
        multicohort_datasets = [ds.name for ds in get_multicohort().get_datasets()]
        sg_ids = []
        sg_vcfs = {}
        for ds in multicohort_datasets:
            sgs = query_for_lrs_vcfs(dataset_name=ds)
            sg_ids.extend(sgs['sg_ids'])
            sg_vcfs.update(sgs['vcfs'])

        assert not set(get_multicohort().get_sequencing_group_ids()) - set(sg_ids), \
            ('SGs in the multicohort do not have VCFs matching the filter criteria: '
             f'{set(get_multicohort().get_sequencing_group_ids()) - set(sg_ids)}'
             ' Adjust the query filter or the input cohorts')

        joint_called = sg_vcfs[sg.id]['meta'].get('joint_called', False)

        # Input VCF and reheadering file
        vcf_path: str = sg_vcfs[sg.id]['vcf']
        lrs_sg_id_mapping = inputs.as_path(get_multicohort(), WriteLrsIdToSgIdMappingFile, 'lrs_sg_id_mapping')

        reformatting_job = bcftools_reformat(
            vcf_path=vcf_path,
            job_name=f'Reformat SNPs Indels VCF for {sg.id}: {"joint-called " if joint_called else ""}{vcf_path}',
            job_attrs={'tool': 'bcftools'},
            lrs_sg_id_mapping_path=lrs_sg_id_mapping,
        )

        outputs = self.expected_outputs(sg)
        get_batch().write_output(reformatting_job.vcf_out, str(outputs['vcf']).removesuffix('.vcf.gz'))

        return self.make_outputs(target=sg, jobs=[reformatting_job], data=outputs)


@stage.stage(required_stages=ModifyVcf)
class MergeVcfsWithBcftools(stage.MultiCohortStage):
    """
    Merge the reformatted SNPs Indels VCFs together with bcftools
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'vcf': self.tmp_prefix / 'snps_indels' / 'merged_reformatted.vcf.gz',
            'index': self.tmp_prefix / 'snps_indels' / 'merged_reformatted.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Use bcftools to merge all the VCFs, and then fill in the tags (requires bcftools 1.18+)
        """
        sgs = get_sgs_from_datasets([d.name for d in multicohort.get_datasets()])

        # Get the reformatted VCFs from the previous stage
        reformatted_vcfs = inputs.as_dict_by_target(ModifyVcf)
        reformatted_vcfs = {
            sg_id: vcf for sg_id, vcf in reformatted_vcfs.items() if sg_id in sgs['vcfs']
        }

        if len(reformatted_vcfs) == 1:
            logger.info('Only one VCF found, skipping merge')
            return None

        vcf_paths = [
            str(reformatted_vcfs[sg_id]['vcf'])
            for sg_id in multicohort.get_sequencing_group_ids()
            if sg_id in reformatted_vcfs
        ]

        merge_job = merge_snps_indels_vcf_with_bcftools(
            batch=get_batch(),
            vcf_paths=vcf_paths,
            job_attrs={'tool': 'bcftools'},
        )

        outputs = self.expected_outputs(multicohort)
        get_batch().write_output(merge_job.output, str(outputs['vcf']).removesuffix('.vcf.gz'))

        return self.make_outputs(multicohort, data=outputs, jobs=merge_job)


@stage.stage(required_stages=[ModifyVcf, MergeVcfsWithBcftools])
class SplitVcfIntoSitesOnlyWithGatk(stage.MultiCohortStage):
    """
    Get the site-only VCFs from the merged VCF by splitting the VCF with GATK SelectVariants,
    then extracting the sites-only VCFs for VEP annotation.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict:
        """
        Generate site-only VCFs from the merged VCF.
        """
        return {
            'siteonly': to_path(self.tmp_prefix / 'snps_indels' / 'siteonly.vcf.gz'),
            'siteonly_part_pattern': str(
                self.tmp_prefix / 'snps_indels' / 'siteonly_parts' / 'part{idx}.vcf.gz'
            ),
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Submit jobs.
        """
        outputs = self.expected_outputs(multicohort)

        scatter_count = joint_calling_scatter_count(len(multicohort.get_sequencing_group_ids()))
        out_siteonly_vcf_part_paths = [
            to_path(outputs['siteonly_part_pattern'].format(idx=idx)) for idx in range(scatter_count)
        ]

        intervals_path = None
        if config_retrieve(['workflow', 'intervals_path'], default=None):
            intervals_path = to_path(config_retrieve(['workflow', 'intervals_path']))

        exclude_intervals_path = None
        if config_retrieve(['workflow', 'exclude_intervals_path'], default=None):
            exclude_intervals_path = to_path(config_retrieve(['workflow', 'exclude_intervals_path']))

        if len(inputs.as_dict_by_target(ModifyVcf)) == 1:
            sg = multicohort.get_sequencing_groups()[0]
            merged_vcf_path = inputs.as_path(sg, ModifyVcf, 'vcf')
        else:
            merged_vcf_path = inputs.as_path(multicohort, MergeVcfsWithBcftools, 'vcf')

        vcf_jobs = split_merged_vcf_and_get_sitesonly_vcfs_for_vep(
            b=get_batch(),
            scatter_count=scatter_count,
            merged_vcf_path=merged_vcf_path,
            tmp_bucket=self.tmp_prefix / 'snps_indels' / 'siteonly_tmp',
            out_siteonly_vcf_part_paths=out_siteonly_vcf_part_paths,
            intervals_path=intervals_path,
            exclude_intervals_path=exclude_intervals_path,
            job_attrs=self.get_job_attrs(),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=vcf_jobs)


@stage.stage(required_stages=[SplitVcfIntoSitesOnlyWithGatk])
class VepLongReadAnnotation(stage.MultiCohortStage):
    """
    Run VEP on the long-read site-only VCFs and write out a Hail table.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort):
        """
        Expected to write a hail table.
        """
        return {'ht': self.tmp_prefix / 'snps_indels' / 'vep.ht'}

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Submit jobs.
        """
        outputs = self.expected_outputs(multicohort)
        scatter_count = joint_calling_scatter_count(len(multicohort.get_sequencing_group_ids()))
        input_siteonly_vcf_part_paths = [
            to_path(
                inputs.as_str(
                    stage=SplitVcfIntoSitesOnlyWithGatk,
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
            tmp_prefix=self.tmp_prefix / 'snps_indels' / 'vep_tmp',
            job_attrs=self.get_job_attrs(),
            scatter_count=scatter_count,
        )
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[
    ModifyVcf,
    MergeVcfsWithBcftools,
    VepLongReadAnnotation]
)
class AnnotateCohortMtFromVcfWithHail(stage.MultiCohortStage):
    """
    First step to transform annotated SNPs Indels callset data into a seqr ready format
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        """
        Expected to write a matrix table.
        """
        return {'mt': self.tmp_prefix / 'snps_indels' / 'cohort.mt'}

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        queue job(s) to rearrange the annotations prior to Seqr transformation
        """
        outputs = self.expected_outputs(multicohort)

        if len(inputs.as_dict_by_target(ModifyVcf)) == 1:
            sg = multicohort.get_sequencing_groups()[0]
            vcf_path = inputs.as_path(sg, ModifyVcf, 'vcf')
        else:
            vcf_path = inputs.as_path(target=multicohort, stage=MergeVcfsWithBcftools, key='vcf')

        vep_ht_path = inputs.as_path(target=multicohort, stage=VepLongReadAnnotation, key='ht')

        job = annotate_cohort_jobs_snps_indels(
            vcf_path=vcf_path,
            out_mt_path=outputs['mt'],
            vep_ht_path=vep_ht_path,
            checkpoint_prefix=self.tmp_prefix / 'snps_indels' / 'checkpoints',
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage.stage(
    required_stages=[AnnotateCohortMtFromVcfWithHail],
    analysis_type='matrixtable',
    analysis_keys=['mt']
)
class SubsetMtToDatasetWithHail(stage.DatasetStage):
    """
    Subset the multicohort Matrixtable to a single dataset
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict:
        """
        Expected to generate a matrix table
        """
        sg_hash = workflow.get_workflow().output_version
        mt_name = f'LongReadSNPsIndels-{sg_hash}-{dataset.name}'
        if family_sgs := get_family_sequencing_groups(dataset):
            return {
                'mt': (dataset.prefix() / 'mt' / f'{mt_name}-{family_sgs["name_suffix"]}.mt'),
            }
        return {
            'mt': (dataset.prefix() / 'mt' / f'{mt_name}.mt'),
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
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

        if family_sgs := get_family_sequencing_groups(dataset):
            sg_ids = family_sgs['family_sg_ids']
        else:
            sg_ids = dataset.get_sequencing_group_ids()

        mt_path = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortMtFromVcfWithHail, key='mt')

        outputs = self.expected_outputs(dataset)

        sg_hash = workflow.get_workflow().output_version
        checkpoint_prefix = dataset.tmp_prefix() / sg_hash / 'snps_indels' / 'mt' / 'checkpoints'

        jobs = annotate_dataset_jobs(
            mt_path=mt_path,
            sg_ids=sg_ids,
            out_mt_path=outputs['mt'],
            tmp_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage.stage(
    required_stages=[SubsetMtToDatasetWithHail],
    analysis_type='es-index',
    analysis_keys=['done_flag'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'SNV_INDEL'},  # noqa: ARG005
)
class ExportSnpsIndelsMtToESIndex(stage.DatasetStage):
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
        index_name = f'{dataset.name}-{sequencing_type}-LRS-SNPsIndels-{sg_hash}'.lower()
        if family_sgs := get_family_sequencing_groups(dataset):
            index_name += f'-{family_sgs["name_suffix"]}'
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'snps_indels' / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
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
        mt_path = inputs.as_str(target=dataset, stage=SubsetMtToDatasetWithHail, key='mt')

        # get the expected outputs as Strings
        index_name = outputs['index_name']
        flag_name = str(outputs['done_flag'])
        # and just the name, used after localisation
        mt_name = mt_path.split('/')[-1]

        req_storage = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(get_sgs_from_datasets([dataset.name])),
        )
        # set all job attributes in one bash
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
