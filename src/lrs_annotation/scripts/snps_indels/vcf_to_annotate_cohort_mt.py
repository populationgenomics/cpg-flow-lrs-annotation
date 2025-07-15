"""
Hail Query functions for seqr loader.
"""
from argparse import ArgumentParser
import logging

import hail as hl

from cpg_utils.config import config_retrieve, reference_path
from cpg_utils.hail_batch import genome_build, init_batch
from cpg_flow.utils import checkpoint_hail
from lrs_annotation.hail_scripts.computed_fields import variant_id, vep
from lrs_annotation.utils import get_init_batch_args_for_job


def annotate_cohort(
    vcf_path: str,
    out_mt_path: str,
    vep_ht_path: str,
    checkpoint_prefix: str | None = None,
    remove_invalid_contigs: bool = False,
):
    """
    Convert VCF to matrix table, annotate for Seqr Loader, add VEP and VQSR
    annotations.
    """
    init_batch(**get_init_batch_args_for_job('annotate_cohort_snps_indels'))

    # tune the logger correctly
    logging.getLogger().setLevel(logging.INFO)

    # hail.zulipchat.com/#narrow/stream/223457-Hail-Batch-support/topic/permissions.20issues/near/398711114
    # don't override the block size, as it explodes the number of partitions when processing TB+ datasets
    # Each partition comes with some computational overhead, it's to be seen whether the standard block size
    # is viable for QOB in large datasets... Requires a test
    mt = hl.import_vcf(
        vcf_path,
        reference_genome=genome_build(),
        skip_invalid_loci=True,
        force_bgz=True,
        array_elements_required=False,
    )
    mt.checkpoint(output=str(checkpoint_prefix) + 'mt-imported.mt', overwrite=True)
    logging.info(f'Imported VCF {vcf_path} as {mt.n_partitions()} partitions')

    # Annotate VEP. Do it before splitting multi, because we run VEP on unsplit VCF,
    # and hl.split_multi_hts can handle multiallelic VEP field.
    vep_ht = hl.read_table(vep_ht_path)
    logging.info(
        f'Adding VEP annotations into the Matrix Table from {vep_ht_path}.'
        f'VEP loaded as {vep_ht.n_partitions()} partitions',
    )
    mt = mt.annotate_rows(vep=vep_ht[mt.locus].vep)

    # Remove any contigs not in the 22 autosomes, X, Y, M
    if remove_invalid_contigs:
        logging.info('Removing invalid contigs')
        # fmt: off
        chromosomes = [
            '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
            '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M',
        ]
        # fmt: on
        standard_contigs = hl.literal(
            hl.literal([f'chr{chrom}' for chrom in chromosomes]),
        )
        mt = mt.filter_rows(standard_contigs.contains(mt.locus.contig))

    # Splitting multi-allelics. We do not handle AS info fields here - we handle
    # them when loading VQSR instead, and populate entire "info" from VQSR.
    mt = hl.split_multi_hts(mt.annotate_rows(locus_old=mt.locus, alleles_old=mt.alleles))
    mt = checkpoint_hail(mt, 'mt-vep-split.mt', checkpoint_prefix)

    ref_ht = hl.read_table(reference_path('seqr_combined_reference_data'))
    clinvar_ht = hl.read_table(reference_path('seqr_clinvar'))

    # If the 'AF' field is already in the entries, we should drop it
    if 'AF' in list(mt.entry):
        logging.info('AF field already present in the entries, dropping it')
        mt = mt.drop('AF')
    logging.info('Annotating with seqr-loader fields: round 1')

    logging.info('Annotating with clinvar and munging annotation fields')
    mt = mt.annotate_rows(
        # still taking just a single value here for downstream compatibility in Seqr
        AC=mt.info.AC[0],
        AF=mt.info.AF[0],
        AN=mt.info.AN,
        aIndex=mt.a_index,
        wasSplit=mt.was_split,
        originalAltAlleles=variant_id.get_expr_for_variant_ids(mt.locus_old, mt.alleles_old),
        sortedTranscriptConsequences=vep.get_expr_for_vep_sorted_transcript_consequences_array(mt.vep),
        variantId=variant_id.get_expr_for_variant_id(mt),
        contig=variant_id.get_expr_for_contig(mt.locus),
        pos=mt.locus.position,
        start=mt.locus.position,
        end=mt.locus.position + hl.len(mt.alleles[0]) - 1,
        ref=mt.alleles[0],
        alt=mt.alleles[1],
        xpos=variant_id.get_expr_for_xpos(mt.locus),
        xstart=variant_id.get_expr_for_xpos(mt.locus),
        xstop=variant_id.get_expr_for_xpos(mt.locus) + hl.len(mt.alleles[0]) - 1,
        clinvar_data=clinvar_ht[mt.row_key],
        ref_data=ref_ht[mt.row_key],
    )

    # this was previously executed in the MtToEs job, as it wasn't possible on QoB
    logging.info('Adding GRCh37 coords')
    liftover_path = reference_path('liftover_38_to_37')
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg38.add_liftover(liftover_path, rg37)
    mt = mt.annotate_rows(rg37_locus=hl.liftover(mt.locus, 'GRCh37'))

    # only remove InbreedingCoeff if present (post-VQSR)
    if 'InbreedingCoeff' in mt.info:
        mt = mt.annotate_rows(info=mt.info.drop('InbreedingCoeff'))

    logging.info(
        'Annotating with seqr-loader fields: round 2 '
        '(expanding sortedTranscriptConsequences, ref_data, clinvar_data)',
    )
    mt = mt.annotate_rows(
        domains=vep.get_expr_for_vep_protein_domains_set_from_sorted(mt.sortedTranscriptConsequences),
        transcriptConsequenceTerms=vep.get_expr_for_vep_consequence_terms_set(mt.sortedTranscriptConsequences),
        transcriptIds=vep.get_expr_for_vep_transcript_ids_set(mt.sortedTranscriptConsequences),
        mainTranscript=vep.get_expr_for_worst_transcript_consequence_annotations_struct(
            mt.sortedTranscriptConsequences,
        ),
        geneIds=vep.get_expr_for_vep_gene_ids_set(mt.sortedTranscriptConsequences),
        codingGeneIds=vep.get_expr_for_vep_gene_ids_set(mt.sortedTranscriptConsequences, only_coding_genes=True),
        cadd=mt.ref_data.cadd,
        dbnsfp=mt.ref_data.dbnsfp,
        geno2mp=mt.ref_data.geno2mp,
        gnomad_exomes=mt.ref_data.gnomad_exomes,
        gnomad_exome_coverage=mt.ref_data.gnomad_exome_coverage,
        gnomad_genomes=mt.ref_data.gnomad_genomes,
        gnomad_genome_coverage=mt.ref_data.gnomad_genome_coverage,
        eigen=mt.ref_data.eigen,
        exac=mt.ref_data.exac,
        g1k=mt.ref_data.g1k,
        mpc=mt.ref_data.mpc,
        primate_ai=mt.ref_data.primate_ai,
        splice_ai=mt.ref_data.splice_ai,
        topmed=mt.ref_data.topmed,
        clinvar=hl.struct(
            allele_id=mt.clinvar_data.info.ALLELEID,
            clinical_significance=hl.delimit(mt.clinvar_data.info.CLNSIG),
            gold_stars=mt.clinvar_data.gold_stars,
        ),
    )
    mt = mt.annotate_globals(
        sourceFilePath=vcf_path,
        genomeVersion=genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
    )
    if sequencing_type := config_retrieve(['workflow', 'sequencing_type']):
        # Map to Seqr-style string
        # https://github.com/broadinstitute/seqr/blob/e0c179c36c0f68c892017de5eab2e4c1b9ffdc92/seqr/models.py#L592-L594
        mt = mt.annotate_globals(
            sampleType={
                'genome': 'WGS',
                'exome': 'WES',
                'single_cell': 'RNA',
            }.get(sequencing_type, ''),
        )

    logging.info('Done:')
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    logging.info(f'Written final matrix table into {out_mt_path}')


def cli_main():
    """
    command line entrypoint
    """

    parser = ArgumentParser()
    parser.add_argument('--vcf_path', type=str, required=True, help='Path to the input VCF file')
    parser.add_argument('--out_mt_path', type=str, required=True, help='Path to the output Matrix Table file')
    parser.add_argument('--vep_ht_path', type=str, required=True, help='Path to the VEP Hail Table file')
    parser.add_argument('--checkpoint_prefix', type=str, default=None, help='Prefix for checkpoint files')
    parser.add_argument('--remove_invalid_contigs', action='store_true', default=False,
                        help='Whether to remove invalid contigs from the Matrix Table')

    args = parser.parse_args()

    annotate_cohort(
        vcf_path=args.vcf_path,
        out_mt_path=args.out_mt_path,
        vep_ht_path=args.vep_ht_path,
        checkpoint_prefix=args.checkpoint_prefix,
        remove_invalid_contigs=args.remove_invalid_contigs,
    )


if __name__ == '__main__':
    cli_main()
