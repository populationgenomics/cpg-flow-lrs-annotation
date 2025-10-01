import gzip
from argparse import ArgumentParser

def cli_main():
    parser = ArgumentParser(description='CLI for the SV VCF modification script')
    parser.add_argument('--vcf_path', help='Path to a localised VCF, this will be modified', required=True)
    parser.add_argument('--vcf_out', help='Path to an output location for the modified VCF', required=True)
    args = parser.parse_args()
    modify_vcf(
        file_in=args.vcf_path,
        file_out=args.vcf_out,
    )

def modify_vcf(file_in: str, file_out: str):
    """
    Scrolls through the VCF and replaces soft-masked characters in the REF allele

    rebuilds the VCF following those edits, and writes the compressed data back out

    Args:
        file_in (str): localised VCF
        file_out (str): local batch output path, same VCF with alterations
    """

    # read and write compressed
    with gzip.open(file_in, 'rt') as f, gzip.open(file_out, 'wt') as f_out:

        for index, line in enumerate(f):

            if index % 10000 == 0:
                print(f'Lines processed: {index}')

            if line.startswith('#'):
                # header line, just write it out unchanged
                continue

            # for non-header lines, split on tabs to extract fields:
            # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE
            l_split = line.rstrip().split('\t')
            l_split[3] = l_split[3].upper()

            # rebuild the string and write as output
            f_out.write('\t'.join(l_split) + '\n')


if __name__ == '__main__':
    cli_main()
