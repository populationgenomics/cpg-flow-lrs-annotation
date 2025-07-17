import hail as hl
from argparse import ArgumentParser

from loguru import logger
from cpg_utils.hail_batch import init_batch
from lrs_annotation.utils import get_init_batch_args_for_job

def annotate_dataset_mt(mt_path: str, out_mt_path: str):
    """
    load the stuff specific to samples in this dataset
    do this after subsetting to specific samples

    Removing the current logic around comparing genotypes to a previous
    callset - doesn't fit with the current implementation

    Args:
        mt_path (str): path to the annotated MatrixTable
        out_mt_path (str): and where do you want it to end up?
    """
    init_batch(**get_init_batch_args_for_job('annotate_dataset_sv'))
    logger.info('Annotating genotypes')

    mt = hl.read_matrix_table(mt_path)
    is_called = hl.is_defined(mt.GT)
    num_alt = hl.if_else(is_called, mt.GT.n_alt_alleles(), -1)
    mt = mt.annotate_rows(
        genotypes=hl.agg.collect(
            hl.struct(
                sample_id=mt.s,
                gq=mt.GQ,
                cn=mt.RD_CN,
                num_alt=num_alt,
            ),
        ),
    )

    def _genotype_filter_samples(fn) -> hl.expr.SetExpression:
        """Filter on the genotypes."""
        return hl.set(mt.genotypes.filter(fn).map(lambda g: g.sample_id))

    # top level - decorator
    def _capture_i_decorator(func):  # noqa: ANN202
        """Call the returned_function(i) which locks in the value of i"""
        def _inner_filter(i):  # noqa: ANN202
            """The _genotype_filter_samples will call this _func with g"""
            def _func(g):  # noqa: ANN202
                """The actual filter function that will be called with g"""
                return func(i, g)

            return _func

        return _inner_filter

    @_capture_i_decorator
    def _filter_num_alt(i, g):  # noqa: ANN202
        return i == g.num_alt

    @_capture_i_decorator
    def _filter_samples_gq(i, g):  # noqa: ANN202
        return (g.gq >= i) & (g.gq < i + 10)

    # github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/sv_mt_schema.py#L221
    # github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/seqr_mt_schema.py#L251
    mt = mt.annotate_rows(
        samples_no_call=_genotype_filter_samples(lambda g: g.num_alt == -1),
        samples_num_alt=hl.struct(**{'%i' % i: _genotype_filter_samples(_filter_num_alt(i)) for i in range(1, 3, 1)}),
        samples_gq_sv=hl.struct(
            **{('%i_to_%i' % (i, i + 10)): _genotype_filter_samples(_filter_samples_gq(i)) for i in range(0, 90, 10)},
        ),
    )

    logger.info('Genotype fields annotated')
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    logger.info(f'Wrote SV MT to {out_mt_path}')



def cli_main():
    """
    command line entrypoint
    """

    parser = ArgumentParser()
    parser.add_argument('--mt_path', type=str, required=True, help='Path to the input MatrixTable file')
    parser.add_argument('--out_mt_path', type=str, required=True, help='Path to write the output MatrixTable file')

    args = parser.parse_args()

    annotate_dataset_mt(
        mt_path=args.mt_path,
        out_mt_path=args.out_mt_path,
    )


if __name__ == '__main__':
    cli_main()
