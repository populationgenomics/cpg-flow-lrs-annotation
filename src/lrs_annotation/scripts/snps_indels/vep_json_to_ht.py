"""
Query function to parse JSON VEP results.

Decoding a JSON string is a very expensive operation, so we need to
provide a fixed schema

This schema is completely rigid (each element needs a name and type),
which causes some issues with parsing, e.g. frequencies:

- the frequencies field is a Dict[str, Dict[str, Float]]
- the top level key is the alt allele, which can be A/C/G/T, or
  "-" for a deletion, and any alt allele for a complex var / indel
- it is not possible to define a fixed schema for this, because any
  alt allele is possible - we have to skip this whole section in the schema
- this isn't an issue for us, as we apply frequencies from an alternative
  resource bundle
- if we were relying on the VEP-supplied frequencies we would have to use
  an alternative way of parsing the VEP JSON, or export to a different format

The key here is that we are aware of the fields that can/will fail, and
avoid relying on them downstream
"""
from argparse import ArgumentParser
import hail as hl
from cpg_utils.hail_batch import init_batch
from lrs_annotation.utils import get_init_batch_args_for_job

def vep_json_to_ht(vep_results_paths: list[str], out_path):
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table.
    receives a vep_version str to determine how data is decoded
    """
    init_batch(**get_init_batch_args_for_job('vep_json_to_ht'))

    # Defining schema inside the function, so we can submit
    # the function to the Batch Job:
    json_schema = hl.dtype(
        """
        struct{
            minimised:int32,
            assembly_name:str,
            allele_string:str,
            ancestral:str,
            colocated_variants:array<struct{
                allele_string:str,
                clin_sig:array<str>,
                clin_sig_allele:str,
                end:int32,
                id:str,
                minimised:int32,
                minor_allele:str,
                minor_allele_freq:float64,
                phenotype_or_disease:int32,
                pubmed:array<int32>,
                seq_region_name:str,
                somatic:int32,
                start:int32,
                strand:int32
            }>,
            context:str,
            end:int32,
            id:str,
            input:str,
            intergenic_consequences:array<struct{
                allele_num:int32,
                consequence_terms:array<str>,
                impact:str,minimised:int32,
                variant_allele:str
            }>,
            most_severe_consequence:str,
            motif_feature_consequences:array<struct{
                allele_num:int32,
                consequence_terms:array<str>,
                high_inf_pos:str,
                impact:str,
                minimised:int32,
                motif_feature_id:str,
                motif_name:str,
                motif_pos:int32,
                motif_score_change:float64,
                strand:int32,
                transcription_factors:array<str>,
                variant_allele:str
            }>,
            regulatory_feature_consequences:array<struct{
                allele_num:int32,
                biotype:str,
                consequence_terms:array<str>,
                impact:str,
                minimised:int32,
                regulatory_feature_id:str,
                variant_allele:str
            }>,
            seq_region_name:str,
            start:int32,
            strand:int32,
            transcript_consequences:array<struct{
                allele_num:int32,
                amino_acids:str,
                appris:str,
                biotype:str,
                canonical:int32,
                mane_select:str,
                mane_plus_clinical:str,
                ccds:str,
                cdna_start:int32,
                cdna_end:int32,
                cds_end:int32,
                cds_start:int32,
                codons:str,
                consequence_terms:array<str>,
                distance:int32,
                domains:array<struct{
                    db:str,
                    name:str
                }>,
                exon:str,
                gene_id:str,
                gene_pheno:int32,
                gene_symbol:str,
                gene_symbol_source:str,
                hgnc_id:str,
                hgvsc:str,
                hgvsp:str,
                hgvs_offset:int32,
                impact:str,
                intron:str,
                lof:str,
                lof_flags:str,
                lof_filter:str,
                lof_info:str,
                existing_inframe_oorfs:int32,
                existing_outofframe_oorfs:int32,
                existing_uorfs:int32,
                5utr_consequence:str,
                5utr_annotation:dict<
                    str,
                    struct{
                        type:str,
                        KozakContext:str,
                        KozakStrength:str,
                        DistanceToCDS:str,
                        CapDistanceToStart:str,
                        DistanceToStop:str,
                        Evidence:str,
                        AltStop:str,
                        AltStopDistanceToCDS:str,
                        FrameWithCDS:str,
                        StartDistanceToCDS:str,
                        newSTOPDistanceToCDS:str,
                        alt_type:str,
                        alt_type_length:str,
                        ref_StartDistanceToCDS:str,
                        ref_type:str,
                        ref_type_length:str
                    }
                >,
                minimised:int32,
                mirna:array<str>,
                polyphen_prediction:str,
                polyphen_score:float64,
                protein_end:int32,
                protein_start:int32,
                protein_id:str,
                sift_prediction:str,
                sift_score:float64,
                strand:int32,
                swissprot:array<str>,
                transcript_id:str,
                trembl:array<str>,
                tsl:int32,
                uniparc:array<str>,
                uniprot_isoform:array<str>,
                variant_allele:str,
                am_class:str,
                am_pathogenicity:float64,
                source:str,
                flags:array<str>
            }>,
            variant_class:str
        }
    """,
    )
    ht = hl.import_table(paths=vep_results_paths, no_header=True, types={'f0': json_schema})
    ht = ht.transmute(vep=ht.f0)
    # Can't use ht.vep.start for start because it can be modified by VEP (e.g. it
    # happens for indels). So instead parsing POS from the original VCF line stored
    # as ht.vep.input field.
    original_vcf_line = ht.vep.input
    start = hl.parse_int(original_vcf_line.split('\t')[1])
    chrom = ht.vep.seq_region_name
    ht = ht.annotate(locus=hl.locus(chrom, start))
    ht = ht.key_by(ht.locus)
    ht.write(str(out_path), overwrite=True)


def cli_main():
    """
    command line entrypoint
    """

    parser = ArgumentParser()
    parser.add_argument('--vep_results_paths', type=str, nargs='+',
                        help='Paths to VEP result files in JSON format')
    parser.add_argument('--out_path', type=str, required=True,
                        help='Output path for the Hail Table')

    args = parser.parse_args()

    vep_json_to_ht(
        vep_results_paths=args.vep_results_paths,
        out_path=args.out_path
    )


if __name__ == '__main__':
    cli_main()
