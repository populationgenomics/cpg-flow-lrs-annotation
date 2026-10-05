"""
Tests for the LongTR pathogenic screening logic.

Threshold fixtures below are the real STRchive values for each locus, inlined rather than
read from references/STRchive-loci.json so the tests do not depend on an untracked file.
"""

import shutil
import subprocess

import pytest
from cyvcf2 import VCF, Variant

from lrs_annotation.scripts.longtr_pathogenic import (
    MATCH_TOLERANCE,
    _add_missing_loci,
    _matches_locus,
    _process_vcf_record,
    _summarise_results,
    _threshold_parts,
    classify_allele,
    compute_allele_repeat_units,
    count_motif_in_sequence,
    scan_vcf,
)

# Standard expansion disorder
HTT = {
    'benign_min': 6,
    'benign_max': 26,
    'intermediate_min': 27,
    'intermediate_max': 35,
    'pathogenic_min': 36,
    'pathogenic_max': 250,
    'inheritance': ['AD'],
    'motif_len': 3,
}

# Discrete/contraction loci: pathogenic range sits below the benign range
VWA1 = {
    'benign_min': 2,
    'benign_max': 2,
    'intermediate_min': None,
    'intermediate_max': None,
    'pathogenic_min': 1,
    'pathogenic_max': 3,
    'inheritance': ['AR'],
    'motif_len': 10,
}
MIR7_2 = {
    'benign_min': 4,
    'benign_max': 4,
    'intermediate_min': None,
    'intermediate_max': None,
    'pathogenic_min': 3,
    'pathogenic_max': 3,
    'inheritance': ['AD'],
    'motif_len': 4,
}

# Repeat counts shared between the synthetic VCF records and the assertions about them
HTT_REF_COUNT, HTT_ALT_COUNT = 18, 44
AR_REF_COUNT, AR_ALT_COUNT = 20, 48

# chrX loci: AR/DMD are X-linked recessive, FMR1 is X-linked dominant
AR = {
    'benign_min': 9,
    'benign_max': 34,
    'intermediate_min': 36,
    'intermediate_max': 37,
    'pathogenic_min': 38,
    'pathogenic_max': 68,
    'inheritance': ['XR'],
    'motif_len': 3,
}
FMR1 = {
    'benign_min': 5,
    'benign_max': 44,
    'intermediate_min': 45,
    'intermediate_max': 200,
    'pathogenic_min': 201,
    'pathogenic_max': 2000,
    'inheritance': ['XD'],
    'motif_len': 3,
}


@pytest.mark.parametrize(
    ('count', 'expected'),
    [
        (20, 'normal'),  # below benign_max
        (26, 'normal'),  # on benign_max
        (27, 'intermediate'),  # on intermediate_min
        (35, 'intermediate'),  # on intermediate_max
        (36, 'pathogenic'),  # on pathogenic_min
        (500, 'pathogenic'),  # beyond pathogenic_max
    ],
)
def test_classify_allele_expansion(count, expected):
    """Expansion thresholds are inclusive at every boundary."""
    assert classify_allele(count, HTT) == expected


def test_classify_allele_gap_between_benign_and_pathogenic():
    """A count above benign_max but below pathogenic_min is intermediate, not uncertain."""
    locus = {'benign_max': 10, 'pathogenic_min': 20, 'intermediate_min': None, 'intermediate_max': None}
    assert classify_allele(15, locus) == 'intermediate'


def test_classify_allele_uncertain_without_pathogenic_threshold():
    """Above benign_max with no pathogenic threshold to compare against is uncertain."""
    locus = {'benign_max': 10, 'pathogenic_min': None, 'intermediate_min': None, 'intermediate_max': None}
    assert classify_allele(25, locus) == 'uncertain'


@pytest.mark.parametrize(
    ('locus', 'count', 'expected'),
    [
        # VWA1: benign is exactly 2; pathogenic 1-3 overlaps it, so benign must win at 2
        (VWA1, 1, 'pathogenic'),
        (VWA1, 2, 'normal'),
        (VWA1, 3, 'pathogenic'),
        (VWA1, 7, 'uncertain'),
        # MIR7-2: 3 pathogenic, 4 normal, anything else unknown
        (MIR7_2, 3, 'pathogenic'),
        (MIR7_2, 4, 'normal'),
        (MIR7_2, 5, 'uncertain'),
    ],
)
def test_classify_allele_discrete(locus, count, expected):
    """Contraction loci check benign before pathogenic, since the ranges overlap."""
    assert classify_allele(count, locus) == expected


@pytest.mark.parametrize(
    ('sequence', 'motif', 'expected'),
    [
        ('CAGCAGCAG', 'CAG', 3),
        ('CAGCAGTCAG', 'CAG', 3),  # interrupting base is not counted as a repeat
        ('', 'CAG', 0),
        ('cagcag', 'CAG', 2),  # case insensitive
        ('AAGGCAAGGC', 'AANGC', 2),  # IUPAC N matches any base
    ],
)
def test_count_motif_in_sequence(sequence, motif, expected):
    """Motif counting is non-overlapping, case-insensitive and IUPAC-aware."""
    assert count_motif_in_sequence(sequence, motif) == expected


def test_compute_allele_repeat_units_selects_ref_and_alt():
    """Allele index 0 counts the ref sequence, non-zero indexes into alt_alleles."""
    ref = 'CAG' * 10
    alt = 'CAG' * 40
    assert compute_allele_repeat_units([0, 1], [alt], ref, 'CAG') == (10.0, 40.0)
    assert compute_allele_repeat_units([0, 0], [alt], ref, 'CAG') == (10.0, 10.0)
    # cyvcf2 reports a '.' ALT as an empty list, so an out-of-range index falls back to ref
    assert compute_allele_repeat_units([1, 1], [], ref, 'CAG') == (10.0, 10.0)


def test_matches_locus_respects_tolerance():
    """A variant span matches a BED entry within MATCH_TOLERANCE bp at both ends, not beyond.

    Offsets are derived from the constant so the tolerance can be retuned without
    editing the expectations.
    """
    start, end = 1000, 1050
    entry = _bed_match('chr1', start, 'CAG', 'L1')
    tol = MATCH_TOLERANCE

    assert _matches_locus(entry, start, end) is True
    assert _matches_locus(entry, start + tol, end + tol) is True  # exactly at tolerance, both ends
    assert _matches_locus(entry, start - tol, end - tol) is True  # and in the other direction
    assert _matches_locus(entry, start + tol + 1, end) is False  # one past, at the start
    assert _matches_locus(entry, start, end + tol + 1) is False  # one past, at the end


VCF_HEADER = [
    '##fileformat=VCFv4.1',
    '##INFO=<ID=START,Number=1,Type=Integer,Description="Repeat start">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="Repeat end">',
    '##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Motif length">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##contig=<ID=chr1,length=250000000>',
    '##contig=<ID=chr4,length=250000000>',
    '##contig=<ID=chrX,length=160000000>',
]


def _vcf_text(rows: list[str], sample: str = 'TESTSAMPLE') -> str:
    """Assemble a valid single-sample VCF from pre-built data rows."""
    chrom_line = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample
    return '\n'.join([*VCF_HEADER, chrom_line, *rows]) + '\n'


def _vcf_row(chrom, start, motif, *, n_ref, n_alt, gt) -> str:
    """Build one VCF data row for a repeat locus spanning 50bp from start."""
    info = f'START={start};END={start + 50};PERIOD={len(motif)}'
    return '\t'.join(
        [chrom, str(start), '.', motif * n_ref, motif * n_alt, '.', 'PASS', info, 'GT:DP', f'{gt}:30'],
    )


def _variant(tmp_path, chrom, start, motif, *, n_ref, n_alt, gt) -> Variant:
    """Write a one-record plain VCF and return it parsed as a cyvcf2 variant.

    Plain text is enough here: _process_vcf_record only needs a variant object, and
    sequential iteration does not require a tabix index.
    """
    path = tmp_path / 'one.vcf'
    path.write_text(_vcf_text([_vcf_row(chrom, start, motif, n_ref=n_ref, n_alt=n_alt, gt=gt)]))
    return next(iter(VCF(str(path))))


def _bed_match(chrom, start, motif, locus_id) -> dict:
    """Build a BED catalog entry for a locus spanning 50bp from start."""
    return {'chrom': chrom, 'start': start, 'end': start + 50, 'motifs': [motif], 'locus_id': locus_id}


@pytest.mark.parametrize(
    ('gt', 'expected_indices', 'expected_gt'),
    [
        ('0|1', [0, 1], '0|1'),
        ('1/1', [1, 1], '1/1'),
        # cyvcf2 puts the phase flag last, so a haploid call must not read it as a second allele.
        # LongTR emits these once --haploid-chrs is supplied for male sex chromosomes.
        ('1', [1, 1], '1'),
        ('0', [0, 0], '0'),
        ('.', [0, 0], '.'),
    ],
)
def test_genotype_parsing_handles_haploid_and_missing(tmp_path, gt, expected_indices, expected_gt):
    """Haploid, diploid and missing calls all resolve to a usable allele pair."""
    variant = _variant(tmp_path, 'chrX', 67545316, 'CAG', n_ref=20, n_alt=28, gt=gt)
    match = _bed_match('chrX', 67545316, 'CAG', 'SBMA_AR')
    result = _process_vcf_record(variant, match, AR, 67545316, 67545366, sex='female')

    counts = {0: 20.0, 1: 28.0}
    assert (result['allele1_ru'], result['allele2_ru']) == tuple(counts[i] for i in expected_indices)
    assert result['gt'] == expected_gt


@pytest.mark.parametrize(('a_ref', 'a_alt'), [(20, 28), (28, 20)])
def test_hemizygous_collapse_keeps_larger_count(tmp_path, a_ref, a_alt):
    """Males get one chrX allele; the larger of LongTR's diploid pair is kept either way round."""
    variant = _variant(tmp_path, 'chrX', 67545316, 'CAG', n_ref=a_ref, n_alt=a_alt, gt='0|1')
    match = _bed_match('chrX', 67545316, 'CAG', 'SBMA_AR')
    result = _process_vcf_record(variant, match, AR, 67545316, 67545366, sex='male')

    assert result['hemizygous'] is True
    assert result['allele1_ru'] == max(a_ref, a_alt)
    assert result['allele2_ru'] is None
    assert result['allele2_status'] is None
    assert result['allele2_seq'] is None


def test_hemizygous_not_applied_to_autosomes(tmp_path):
    """An autosomal locus keeps both alleles even for a male sample."""
    variant = _variant(tmp_path, 'chr4', 3074876, 'CAG', n_ref=HTT_REF_COUNT, n_alt=HTT_ALT_COUNT, gt='0|1')
    match = _bed_match('chr4', 3074876, 'CAG', 'HD_HTT')
    result = _process_vcf_record(variant, match, HTT, 3074876, 3074926, sex='male')

    assert result['hemizygous'] is False
    assert (result['allele1_ru'], result['allele2_ru']) == (HTT_REF_COUNT, HTT_ALT_COUNT)
    assert result['locus_status'] == 'pathogenic'


@pytest.mark.parametrize(
    ('sex', 'expected_status'),
    [
        ('female', 'carrier'),  # one pathogenic allele on an X-linked recessive locus
        ('male', 'pathogenic'),  # hemizygous, so the expansion is expressed
        ('unknown', 'pathogenic'),  # never downgrade without knowing sex
    ],
)
def test_xlinked_recessive_single_pathogenic_allele(tmp_path, sex, expected_status):
    """A female heterozygote at an XR locus is a carrier; males and unknown sex are not."""
    variant = _variant(tmp_path, 'chrX', 67545316, 'CAG', n_ref=AR_REF_COUNT, n_alt=AR_ALT_COUNT, gt='0|1')
    match = _bed_match('chrX', 67545316, 'CAG', 'SBMA_AR')
    result = _process_vcf_record(variant, match, AR, 67545316, 67545366, sex=sex)
    assert result['locus_status'] == expected_status


def test_female_homozygous_xr_is_affected_not_carrier(tmp_path):
    """Two pathogenic alleles at an XR locus means affected, so carrier must not apply."""
    variant = _variant(tmp_path, 'chrX', 67545316, 'CAG', n_ref=AR_REF_COUNT, n_alt=AR_ALT_COUNT, gt='1|1')
    match = _bed_match('chrX', 67545316, 'CAG', 'SBMA_AR')
    result = _process_vcf_record(variant, match, AR, 67545316, 67545366, sex='female')

    assert result['allele1_status'] == result['allele2_status'] == 'pathogenic'
    assert result['locus_status'] == 'pathogenic'


def test_female_xlinked_dominant_stays_pathogenic(tmp_path):
    """X-linked dominant loci manifest in females, so carrier must not apply to FMR1."""
    variant = _variant(tmp_path, 'chrX', 147912049, 'CGG', n_ref=30, n_alt=250, gt='0|1')
    match = _bed_match('chrX', 147912049, 'CGG', 'FXS_FMR1')
    result = _process_vcf_record(variant, match, FMR1, 147912049, 147912099, sex='female')
    assert result['locus_status'] == 'pathogenic'


def test_threshold_parts_enumerates_discrete_loci():
    """Discrete loci list their qualifying counts; expansion loci keep range notation."""
    assert _threshold_parts(VWA1) == [
        'Normal: 2',
        'Pathogenic: 1, 3',
        'any other count: uncertain',
    ]
    assert _threshold_parts(MIR7_2) == [
        'Normal: 4',
        'Pathogenic: 3',
        'any other count: uncertain',
    ]
    assert _threshold_parts(HTT) == [
        'Normal: ≤26',
        'Intermediate: 27-35',
        'Pathogenic: ≥36',
    ]


def test_add_missing_loci_fills_every_field_the_template_reads():
    """Loci absent from the VCF still need the full key set, since the template reads them unconditionally."""
    entries = [
        _bed_match('chr4', 3074876, 'CAG', 'HD_HTT'),
        _bed_match('chrX', 67545316, 'CAG', 'SBMA_AR'),
    ]
    results = {}
    _add_missing_loci(results, entries, {'HD_HTT': HTT, 'SBMA_AR': AR})

    assert set(results) == {'HD_HTT', 'SBMA_AR'}
    entry = results['HD_HTT']
    assert entry['genotyped'] is False
    assert entry['locus_status'] == 'not_genotyped'
    assert entry['allele1_status'] == entry['allele2_status'] == 'not_genotyped'
    assert entry['allele1_ru'] is entry['allele2_ru'] is None
    # hemizygous and thresholds are read for every card, genotyped or not
    assert entry['hemizygous'] is False
    assert entry['thresholds'] == ['Normal: ≤26', 'Intermediate: 27-35', 'Pathogenic: ≥36']


def test_add_missing_loci_leaves_genotyped_entries_alone():
    """An already-genotyped locus must not be overwritten by a not-genotyped placeholder."""
    entries = [_bed_match('chr4', 3074876, 'CAG', 'HD_HTT')]
    results = {'HD_HTT': {'locus_status': 'pathogenic', 'genotyped': True}}
    _add_missing_loci(results, entries, {'HD_HTT': HTT})

    assert results['HD_HTT']['locus_status'] == 'pathogenic'
    assert results['HD_HTT']['genotyped'] is True


def test_summarise_results_counts_by_status():
    """Totals drive the summary tiles and filter button counts."""
    rows = [
        {'genotyped': True, 'locus_status': 'pathogenic'},
        {'genotyped': True, 'locus_status': 'carrier'},
        {'genotyped': True, 'locus_status': 'normal'},
        {'genotyped': True, 'locus_status': 'normal'},
        {'genotyped': False, 'locus_status': 'not_genotyped'},
    ]
    assert _summarise_results(rows) == {
        'total_loci': 5,
        'genotyped': 4,
        'not_genotyped': 1,
        'pathogenic': 1,
        'carrier': 1,
        'normal': 2,
    }


def test_summarise_results_handles_no_results():
    """An empty loci list must not KeyError when the report is rendered."""
    assert _summarise_results([]) == {'total_loci': 0, 'genotyped': 0, 'not_genotyped': 0}


# scan_vcf fetches by region, which htslib can only do against a bgzipped+indexed file
needs_htslib = pytest.mark.skipif(
    shutil.which('bgzip') is None or shutil.which('tabix') is None,
    reason='bgzip/tabix needed to build an indexed VCF for region queries',
)


@pytest.fixture
def indexed_vcf(tmp_path):
    """Write a two-record VCF covering one autosomal and one chrX locus, bgzipped and indexed."""
    bgzip, tabix = shutil.which('bgzip'), shutil.which('tabix')
    plain = tmp_path / 'mini.vcf'
    plain.write_text(
        _vcf_text(
            [
                _vcf_row('chr4', 3074876, 'CAG', n_ref=HTT_REF_COUNT, n_alt=HTT_ALT_COUNT, gt='0|1'),
                _vcf_row('chrX', 67545316, 'CAG', n_ref=AR_REF_COUNT, n_alt=AR_ALT_COUNT, gt='0|1'),
            ],
        ),
    )
    subprocess.run([bgzip, '-f', str(plain)], check=True)  # noqa: S603
    gz = tmp_path / 'mini.vcf.gz'
    subprocess.run([tabix, '-f', '-p', 'vcf', str(gz)], check=True)  # noqa: S603
    return str(gz)


@pytest.fixture
def bed_entries():
    """The flat BED catalog scan_vcf iterates, one entry per disease locus."""
    return [
        _bed_match('chr4', 3074876, 'CAG', 'HD_HTT'),
        _bed_match('chrX', 67545316, 'CAG', 'SBMA_AR'),
    ]


def _write_indexed(tmp_path, rows: list[str]) -> str:
    """bgzip + tabix a VCF built from the given rows, returning the .gz path."""
    bgzip, tabix = shutil.which('bgzip'), shutil.which('tabix')
    plain = tmp_path / 'offset.vcf'
    plain.write_text(_vcf_text(rows))
    subprocess.run([bgzip, '-f', str(plain)], check=True)  # noqa: S603
    gz = tmp_path / 'offset.vcf.gz'
    subprocess.run([tabix, '-f', '-p', 'vcf', str(gz)], check=True)  # noqa: S603
    return str(gz)


@needs_htslib
@pytest.mark.parametrize('n_ref', [3, HTT_REF_COUNT])
def test_scan_vcf_finds_short_ref_offset_behind_the_catalog_entry(tmp_path, n_ref):
    """A record offset backwards by the full tolerance must still be fetched.

    Over a third of the genotyped loci have a reference repeat shorter than the tolerance
    (FXN and PABPN1 are 19bp), so such a record's span can end before the catalog entry
    begins. The region fetch is padded by the tolerance precisely to still reach it.
    """
    entry = _bed_match('chr4', 3074876, 'CAG', 'HD_HTT')
    vcf = _write_indexed(
        tmp_path,
        [_vcf_row('chr4', 3074876 - MATCH_TOLERANCE, 'CAG', n_ref=n_ref, n_alt=n_ref + 1, gt='0|1')],
    )
    results = scan_vcf(vcf, [entry], {'HD_HTT': HTT}, 'female')
    assert results[0]['genotyped'] is True


@needs_htslib
@pytest.mark.parametrize(
    ('offset', 'should_match'),
    [
        (0, True),  # exact
        (-MATCH_TOLERANCE, True),  # at tolerance, record starts before the catalog entry
        (MATCH_TOLERANCE, True),  # at tolerance, record starts after
        (-(MATCH_TOLERANCE + 1), False),  # one past tolerance
        (MATCH_TOLERANCE + 1, False),  # one past tolerance
    ],
)
def test_scan_vcf_tolerates_offset_coordinates(tmp_path, offset, should_match):
    """LongTR coordinates do not exactly match STRchive's, so matching must survive a shift."""
    entry = _bed_match('chr4', 3074876, 'CAG', 'HD_HTT')
    vcf = _write_indexed(
        tmp_path,
        [_vcf_row('chr4', 3074876 + offset, 'CAG', n_ref=HTT_REF_COUNT, n_alt=HTT_ALT_COUNT, gt='0|1')],
    )
    results = scan_vcf(vcf, [entry], {'HD_HTT': HTT}, 'female')
    assert results[0]['genotyped'] is should_match


@needs_htslib
def test_scan_vcf_end_to_end(indexed_vcf, bed_entries):
    """Reads a VCF, matches both loci, and sorts worst-first."""
    strchive = {'HD_HTT': HTT, 'SBMA_AR': AR}
    results = scan_vcf(indexed_vcf, bed_entries, strchive, 'female')

    by_id = {r['locus_id']: r for r in results}
    assert set(by_id) == {'HD_HTT', 'SBMA_AR'}

    assert by_id['HD_HTT']['locus_status'] == 'pathogenic'
    assert (by_id['HD_HTT']['allele1_ru'], by_id['HD_HTT']['allele2_ru']) == (HTT_REF_COUNT, HTT_ALT_COUNT)
    # female + XR + one pathogenic allele
    assert by_id['SBMA_AR']['locus_status'] == 'carrier'

    # pathogenic sorts ahead of carrier
    assert [r['locus_id'] for r in results] == ['HD_HTT', 'SBMA_AR']


@needs_htslib
def test_scan_vcf_male_collapses_chrx(indexed_vcf, bed_entries):
    """The same VCF read as male yields a hemizygous chrX call and an untouched autosome."""
    strchive = {'HD_HTT': HTT, 'SBMA_AR': AR}
    results = scan_vcf(indexed_vcf, bed_entries, strchive, 'male')
    by_id = {r['locus_id']: r for r in results}

    assert by_id['SBMA_AR']['hemizygous'] is True
    assert by_id['SBMA_AR']['allele1_ru'] == AR_ALT_COUNT
    assert by_id['SBMA_AR']['allele2_ru'] is None
    assert by_id['HD_HTT']['hemizygous'] is False
    assert by_id['HD_HTT']['allele2_ru'] == HTT_ALT_COUNT
