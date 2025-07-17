import hail as hl
from argparse import ArgumentParser

from loguru import logger
from cpg_utils.hail_batch import init_batch
from lrs_annotation.utils import get_init_batch_args_for_job

def subset_mt_to_sgs(mt_path: str, sg_ids: list[str], out_mt_path: str):
    """
    Subset the MatrixTable to the provided list of sgs and to variants present
    in those sgs

    Args:
        mt_path (str): cohort-level matrix table from VCF.
        sg_ids (list[str]): sgs to take from the matrix table.
        out_mt_path (str): path to write the result.
    """
    init_batch(**get_init_batch_args_for_job('subset_mt_to_sgs'))
    mt = hl.read_matrix_table(mt_path)

    unique_sg_ids: set[str] = set(sg_ids)

    mt_sg_ids = set(mt.s.collect())

    if sg_ids_not_in_mt := unique_sg_ids - mt_sg_ids:
        raise ValueError(
            f'Found {len(sg_ids_not_in_mt)}/{len(unique_sg_ids)} IDs in the requested subset not in the callset.\n'
            f'IDs that aren\'t in the callset: {sg_ids_not_in_mt}\n'
            f'All callset SG IDs: {mt_sg_ids}',
        )

    logger.info(f'Found {len(mt_sg_ids)} samples in mt, subsetting to {len(unique_sg_ids)} samples.')

    mt = mt.filter_cols(hl.literal(unique_sg_ids).contains(mt.s))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt.write(out_mt_path, overwrite=True)
    logger.info(f'Finished subsetting to {len(unique_sg_ids)} samples, written to {out_mt_path}.')

def cli_main():
    """
    command line entrypoint
    """

    parser = ArgumentParser()
    parser.add_argument('--mt_path', type=str, required=True, help='Path to the input MatrixTable file')
    parser.add_argument('--sg_ids', type=str, required=True, help='List of sequencing group IDs to subset')
    parser.add_argument('--out_mt_path', type=str, required=True, help='Path to write the output MatrixTable file')

    args = parser.parse_args()

    subset_mt_to_sgs(
        mt_path=args.mt_path,
        sg_ids=args.sg_ids.split(','),
        out_mt_path=args.out_mt_path,
    )


if __name__ == '__main__':
    cli_main()
