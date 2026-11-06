"""
Tests for the SV VCF merge job.
"""

import shutil
import subprocess

from lrs_annotation.jobs.svs.MergeVcfs import UNIQUE_ID_AWK

AWK = shutil.which('awk')
VCF_HEADER = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'


def _run_awk(vcf_text: str) -> str:
    assert AWK
    return subprocess.run([AWK, UNIQUE_ID_AWK], input=vcf_text, capture_output=True, text=True, check=True).stdout  # noqa: S603


def _uniquify(body: str) -> list[list[str]]:
    return [line.split('\t') for line in _run_awk(VCF_HEADER + body).splitlines() if not line.startswith('#')]


def test_header_is_unchanged():
    assert _run_awk(VCF_HEADER) == VCF_HEADER


def test_duplicate_ids_are_made_unique_keeping_first():
    body = (
        'chr1\t100\tINS_chr1_100_100\tA\t<INS>\t.\tPASS\tSVLEN=50\n'
        'chr1\t100\tINS_chr1_100_100\tA\t<INS>\t.\tPASS\tSVLEN=80\n'
        'chr1\t100\tINS_chr1_100_100\tA\t<INS>\t.\tPASS\tSVLEN=90\n'
        'chr1\t200\tDEL_chr1_200_300\tC\t<DEL>\t.\tPASS\tEND=300\n'
    )
    records = _uniquify(body)
    assert [r[2] for r in records] == [
        'INS_chr1_100_100',
        'INS_chr1_100_100_2',
        'INS_chr1_100_100_3',
        'DEL_chr1_200_300',
    ]


def test_other_columns_are_preserved():
    body = 'chr1\t100\tX\tA\t<INS>\t.\tPASS\tSVLEN=50;NOTE=a b\nchr1\t100\tX\tA\t<INS>\t.\tPASS\tSVLEN=80\n'
    records = _uniquify(body)
    assert records[0] == ['chr1', '100', 'X', 'A', '<INS>', '.', 'PASS', 'SVLEN=50;NOTE=a b']
    assert records[1] == ['chr1', '100', 'X_2', 'A', '<INS>', '.', 'PASS', 'SVLEN=80']
