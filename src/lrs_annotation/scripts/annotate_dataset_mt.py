import hail as hl
from argparse import ArgumentParser

from loguru import logger


def annotate_dataset_mt(mt_path: str, out_mt_path: str):
    """
    Add dataset-level annotations.
    """
    mt = hl.read_matrix_table(mt_path)

    # Convert the mt genotype entries into num_alt, gq, ab, dp, and sample_id.
    is_called = hl.is_defined(mt.GT)
    genotype_fields = {
        'num_alt': hl.if_else(is_called, mt.GT.n_alt_alleles(), -1),
        'gq': hl.if_else(is_called, mt.GQ, hl.null(hl.tint)),
        'ab': hl.bind(
            lambda total: hl.if_else(
                is_called & (total != 0) & (hl.len(mt.AD) > 1),
                hl.float(mt.AD[1] / total),
                hl.missing(hl.tfloat),
            ),
            hl.sum(mt.AD),
        ),
        'dp': hl.if_else(is_called, hl.int(hl.min(mt.DP, 32000)), hl.missing(hl.tfloat)),
        'sample_id': mt.s,
    }
    logger.info('Annotating genotypes')
    mt = mt.annotate_rows(genotypes=hl.agg.collect(hl.struct(**genotype_fields)))

    def _genotype_filter_samples(fn):  # noqa: ANN202
        """Filter on the genotypes."""
        return hl.set(mt.genotypes.filter(fn).map(lambda g: g.sample_id))

    # Colocated variants are represented as a semicolon-separated list of rsids.
    # We only want the first one, and we want to truncate it to 512 characters.
    # This is to ensure that the rsid field is not too long for the es export.
    mt = mt.annotate_rows(
        rsid=hl.str(
            hl.if_else(
                hl.str(mt.rsid).contains(';'),
                hl.str(mt.rsid).split(';')[0],
                mt.rsid,
            )[:512],
        ),
    )

    # 2022-07-28 mfranklin: Initially the code looked like:
    #           {**_genotype_filter_samples(lambda g: g.num_alt == i) for i in ...}
    #   except the lambda definition doesn't bind the loop variable i in this scope
    #   instead let's define the filters as functions, and wrap them in a decorator
    #   that captures the value of i.

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
        return (g.gq >= i) & (g.gq < i + 5)

    @_capture_i_decorator
    def _filter_samples_ab(i, g):  # noqa: ANN202
        return (g.num_alt == 1) & ((g.ab * 100) >= i) & ((g.ab * 100) < i + 5)

    mt = mt.annotate_rows(
        samples_no_call=_genotype_filter_samples(lambda g: g.num_alt == -1),
        samples_num_alt=hl.struct(**{('%i' % i): _genotype_filter_samples(_filter_num_alt(i)) for i in range(1, 3, 1)}),
        samples_gq=hl.struct(
            **{('%i_to_%i' % (i, i + 5)): _genotype_filter_samples(_filter_samples_gq(i)) for i in range(0, 95, 5)},
        ),
        samples_ab=hl.struct(
            **{'%i_to_%i' % (i, i + 5): _genotype_filter_samples(_filter_samples_ab(i)) for i in range(0, 45, 5)},
        ),
    )
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    logger.info(f'Written {out_mt_path}')


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
