from argparse import ArgumentParser

import hail as hl
from cpg_utils.hail_batch import genome_build, init_batch
from loguru import logger

from lrs_annotation.utils import get_init_batch_args_for_job


def vcf_to_unannotated_mt(
    vcf_path: str,
    out_mt_path: str,
):
    """
    Convert VCF to matrix table, unannotated for Talos.
    """
    init_batch(**get_init_batch_args_for_job('vcf_to_unannotated_mt'))

    mt = hl.import_vcf(
        vcf_path,
        reference_genome=genome_build(),
        skip_invalid_loci=True,
        force_bgz=True,
        array_elements_required=False,
    )
    logger.info(f'Imported VCF {vcf_path} as {mt.n_partitions()} partitions')
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    logger.info(f'Wrote final matrix table to {out_mt_path}')


def cli_main():
    """
    command line entrypoint
    """
    parser = ArgumentParser()
    parser.add_argument('--vcf_path', type=str, required=True, help='Path to the input VCF file')
    parser.add_argument('--out_mt_path', type=str, required=True, help='Path to the output Matrix Table file')
    args = parser.parse_args()

    vcf_to_unannotated_mt(
        vcf_path=args.vcf_path,
        out_mt_path=args.out_mt_path,
    )


if __name__ == '__main__':
    cli_main()
