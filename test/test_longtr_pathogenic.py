"""
Tests for the LongTR pathogenic screening logic.

Threshold fixtures below are the real STRchive values for each locus, inlined rather than
read from references/STRchive-loci.json so the tests do not depend on an untracked file.
"""

import pytest

from lrs_annotation.scripts.longtr_pathogenic import (
    _add_missing_loci,
    _parse_gt_indices,
    _process_vcf_record,
    _summarise_results,
    _threshold_parts,
    classify_allele,
    compute_allele_repeat_units,
    count_motif_in_sequence,
    find_matching_loci,
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


@pytest.mark.parametrize(
    ('gt', 'expected'),
    [
        ('0/1', [0, 1]),
        ('1|2', [1, 2]),
        ('1', [1, 1]),  # haploid call is duplicated - see test_hemizygous_* for why that matters
        ('.', [0, 0]),
        ('./.', [0, 0]),
    ],
)
def test_parse_gt_indices(gt, expected):
    """GT parsing handles phased, unphased, haploid and missing calls."""
    assert _parse_gt_indices(gt) == expected


def test_compute_allele_repeat_units_selects_ref_and_alt():
    """Allele index 0 counts the ref sequence, non-zero indexes into alt_alleles."""
    ref = 'CAG' * 10
    alt = 'CAG' * 40
    assert compute_allele_repeat_units([0, 1], [alt], ref, 'CAG') == (10.0, 40.0)
    assert compute_allele_repeat_units([0, 0], [alt], ref, 'CAG') == (10.0, 10.0)
    # a '.' alt falls back to the ref count rather than erroring
    assert compute_allele_repeat_units([1, 1], ['.'], ref, 'CAG') == (10.0, 10.0)


def test_find_matching_loci_respects_tolerance():
    """BED entries match within MATCH_TOLERANCE bp at both ends, and not beyond it."""
    index = {'chr1': [{'chrom': 'chr1', 'start': 1000, 'end': 1050, 'motifs': ['CAG'], 'locus_id': 'L1'}]}
    assert len(find_matching_loci('chr1', 1000, 1050, index)) == 1
    assert len(find_matching_loci('chr1', 1020, 1070, index)) == 1  # exactly 20 off at both ends
    assert len(find_matching_loci('chr1', 1021, 1050, index)) == 0  # 21 off at the start
    assert len(find_matching_loci('chr1', 1000, 1050, {'chr2': index['chr1']})) == 0  # wrong chromosome


def _vcf_record(chrom, start, motif, *, n_ref, n_alt, gt) -> list[str]:
    """Build the minimal VCF column list that _process_vcf_record consumes."""
    return [chrom, str(start), '.', motif * n_ref, motif * n_alt, '.', 'PASS', '.', 'GT:DP', f'{gt}:30']


def _bed_match(chrom, start, motif, locus_id) -> dict:
    """Build a BED catalog entry for a locus spanning 50bp from start."""
    return {'chrom': chrom, 'start': start, 'end': start + 50, 'motifs': [motif], 'locus_id': locus_id}


@pytest.mark.parametrize(('a_ref', 'a_alt'), [(20, 28), (28, 20)])
def test_hemizygous_collapse_keeps_larger_count(a_ref, a_alt):
    """Males get one chrX allele; the larger of LongTR's diploid pair is kept either way round."""
    cols = _vcf_record('chrX', 67545316, 'CAG', n_ref=a_ref, n_alt=a_alt, gt='0|1')
    match = _bed_match('chrX', 67545316, 'CAG', 'SBMA_AR')
    result = _process_vcf_record(cols, match, AR, 67545316, 67545366, {'PERIOD': '3'}, 'male')

    assert result['hemizygous'] is True
    assert result['allele1_ru'] == max(a_ref, a_alt)
    assert result['allele2_ru'] is None
    assert result['allele2_status'] is None
    assert result['allele2_seq'] is None


def test_hemizygous_not_applied_to_autosomes():
    """An autosomal locus keeps both alleles even for a male sample."""
    cols = _vcf_record('chr4', 3074876, 'CAG', n_ref=HTT_REF_COUNT, n_alt=HTT_ALT_COUNT, gt='0|1')
    match = _bed_match('chr4', 3074876, 'CAG', 'HD_HTT')
    result = _process_vcf_record(cols, match, HTT, 3074876, 3074926, {'PERIOD': '3'}, 'male')

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
def test_xlinked_recessive_single_pathogenic_allele(sex, expected_status):
    """A female heterozygote at an XR locus is a carrier; males and unknown sex are not."""
    cols = _vcf_record('chrX', 67545316, 'CAG', n_ref=AR_REF_COUNT, n_alt=AR_ALT_COUNT, gt='0|1')
    match = _bed_match('chrX', 67545316, 'CAG', 'SBMA_AR')
    result = _process_vcf_record(cols, match, AR, 67545316, 67545366, {'PERIOD': '3'}, sex)
    assert result['locus_status'] == expected_status


def test_female_homozygous_xr_is_affected_not_carrier():
    """Two pathogenic alleles at an XR locus means affected, so carrier must not apply."""
    cols = _vcf_record('chrX', 67545316, 'CAG', n_ref=AR_REF_COUNT, n_alt=AR_ALT_COUNT, gt='1|1')
    match = _bed_match('chrX', 67545316, 'CAG', 'SBMA_AR')
    result = _process_vcf_record(cols, match, AR, 67545316, 67545366, {'PERIOD': '3'}, 'female')

    assert result['allele1_status'] == result['allele2_status'] == 'pathogenic'
    assert result['locus_status'] == 'pathogenic'


def test_female_xlinked_dominant_stays_pathogenic():
    """X-linked dominant loci manifest in females, so carrier must not apply to FMR1."""
    cols = _vcf_record('chrX', 147912049, 'CGG', n_ref=30, n_alt=250, gt='0|1')
    match = _bed_match('chrX', 147912049, 'CGG', 'FXS_FMR1')
    result = _process_vcf_record(cols, match, FMR1, 147912049, 147912099, {'PERIOD': '3'}, 'female')
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
    index = {
        'chr4': [_bed_match('chr4', 3074876, 'CAG', 'HD_HTT')],
        'chrX': [_bed_match('chrX', 67545316, 'CAG', 'SBMA_AR')],
    }
    results = {}
    _add_missing_loci(results, index, {'HD_HTT': HTT, 'SBMA_AR': AR})

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
    index = {'chr4': [_bed_match('chr4', 3074876, 'CAG', 'HD_HTT')]}
    results = {'HD_HTT': {'locus_status': 'pathogenic', 'genotyped': True}}
    _add_missing_loci(results, index, {'HD_HTT': HTT})

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


MINIMAL_VCF = """\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTESTSAMPLE
chr4\t3074876\t.\t{htt_ref}\t{htt_alt}\t.\tPASS\tSTART=3074876;END=3074926;PERIOD=3\tGT:DP\t0|1:30
chrX\t67545316\t.\t{ar_ref}\t{ar_alt}\t.\tPASS\tSTART=67545316;END=67545366;PERIOD=3\tGT:DP\t0|1:30
"""


@pytest.fixture
def minimal_vcf(tmp_path):
    """Write a two-record VCF covering one autosomal and one chrX locus."""
    path = tmp_path / 'mini.vcf'
    path.write_text(
        MINIMAL_VCF.format(
            htt_ref='CAG' * HTT_REF_COUNT,
            htt_alt='CAG' * HTT_ALT_COUNT,
            ar_ref='CAG' * AR_REF_COUNT,
            ar_alt='CAG' * AR_ALT_COUNT,
        )
    )
    return str(path)


@pytest.fixture
def minimal_index():
    return {
        'chr4': [_bed_match('chr4', 3074876, 'CAG', 'HD_HTT')],
        'chrX': [_bed_match('chrX', 67545316, 'CAG', 'SBMA_AR')],
    }


def test_scan_vcf_end_to_end(minimal_vcf, minimal_index):
    """Reads a VCF, matches both loci, and sorts worst-first."""
    strchive = {'HD_HTT': HTT, 'SBMA_AR': AR}
    results = scan_vcf(minimal_vcf, minimal_index, strchive, 'female')

    by_id = {r['locus_id']: r for r in results}
    assert set(by_id) == {'HD_HTT', 'SBMA_AR'}

    assert by_id['HD_HTT']['locus_status'] == 'pathogenic'
    assert (by_id['HD_HTT']['allele1_ru'], by_id['HD_HTT']['allele2_ru']) == (HTT_REF_COUNT, HTT_ALT_COUNT)
    # female + XR + one pathogenic allele
    assert by_id['SBMA_AR']['locus_status'] == 'carrier'

    # pathogenic sorts ahead of carrier
    assert [r['locus_id'] for r in results] == ['HD_HTT', 'SBMA_AR']


def test_scan_vcf_male_collapses_chrx(minimal_vcf, minimal_index):
    """The same VCF read as male yields a hemizygous chrX call and an untouched autosome."""
    strchive = {'HD_HTT': HTT, 'SBMA_AR': AR}
    results = scan_vcf(minimal_vcf, minimal_index, strchive, 'male')
    by_id = {r['locus_id']: r for r in results}

    assert by_id['SBMA_AR']['hemizygous'] is True
    assert by_id['SBMA_AR']['allele1_ru'] == AR_ALT_COUNT
    assert by_id['SBMA_AR']['allele2_ru'] is None
    assert by_id['HD_HTT']['hemizygous'] is False
    assert by_id['HD_HTT']['allele2_ru'] == HTT_ALT_COUNT
