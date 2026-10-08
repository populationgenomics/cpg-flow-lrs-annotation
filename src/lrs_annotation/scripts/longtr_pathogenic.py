"""
Screen a LongTR VCF against STRchive disease-associated TR loci and generate
a standalone HTML report with gauge visualizations.

Algorithm:
1. Load STRchive locus metadata (JSON) and genomic coordinates (BED).
2. Index BED entries by chromosome for fast positional lookup.
3. Scan the VCF: for each variant, find all matching disease loci within a
   positional tolerance and compute per-allele repeat unit counts using
   motif counting for both ref and alt alleles. This handles interruptions
   (bases between motif repeats) correctly, whereas a length-based approach
   would overcount by including interrupting bases as part of the repeat.

4. Classify each allele as normal/intermediate/pathogenic/uncertain against
   STRchive thresholds, then derive a per-locus status from the worst allele.
5. Append not-genotyped entries for any disease loci absent from the VCF.
6. Render results as an interactive HTML report (SVG gauges, motif
   highlighting, evidence filtering) and a structured JSON summary.

References:
  - STRchive-loci.json: disease locus metadata + thresholds
  - STRchive-disease-loci.hg38.longTR.bed: coordinates for matching
  Both from STRchive (github.com/dashnowlab/STRchive), with custom entries.
"""

import datetime
import json
import math
import re
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path

import jinja2
from cyvcf2 import VCF
from loguru import logger
from markupsafe import Markup

MATCH_TOLERANCE = 20

# Worst-first ordering, used both to reduce two alleles to a locus call and to sort the report.
STATUS_PRIORITY = {
    'pathogenic': 0,
    'intermediate': 1,
    'uncertain': 2,
    'normal': 3,
    'not_genotyped': 4,
}

INHERITANCE_NAMES = {
    'AD': 'Autosomal dominant',
    'AR': 'Autosomal recessive',
    'XD': 'X-linked dominant',
    'XR': 'X-linked recessive',
}

EXTERNAL_LINK_DEFS = [
    ('omim', 'OMIM', 'https://omim.org/entry/{}'),
    ('genereviews', 'GeneReviews', 'https://www.ncbi.nlm.nih.gov/books/{}'),
    ('gnomad', 'gnomAD', 'https://gnomad.broadinstitute.org/short-tandem-repeat/{}?dataset=gnomad_r4'),
    ('stripy', 'STRipy', 'https://stripy.org/database/{}'),
    ('medgen', 'MedGen', 'https://www.ncbi.nlm.nih.gov/medgen/{}'),
    ('gard', 'GARD', 'https://rarediseases.info.nih.gov/diseases/{}'),
    ('orphanet', 'Orphanet', 'https://www.orpha.net/en/disease/detail/{}'),
]


def load_strchive_json(path: str, custom_path: str | None = None) -> dict:
    """Load the STRchive catalog by locus ID, overlaying locally curated loci on top.

    The base file is kept as a pristine copy of upstream STRchive so it can be re-synced
    wholesale. Anything local - extra loci, or thresholds upstream does not publish - lives in
    the overlay and wins on an ID collision. Overlay entries may be partial, in which case the
    remaining fields are inherited from upstream.
    """
    with open(path) as f:
        loci = {locus['id']: locus for locus in json.load(f)}
    upstream_count = len(loci)

    if custom_path:
        with open(custom_path) as f:
            custom = json.load(f)
        added = sum(1 for entry in custom if entry['id'] not in loci)
        for entry in custom:
            locus_id = entry['id']
            loci[locus_id] = {**loci.get(locus_id, {}), **entry}
        logger.info(
            f'Loaded {upstream_count} loci from STRchive, overlaid {len(custom)} locally curated '
            f'({added} new, {len(custom) - added} overriding upstream)',
        )

    return loci


def catalog_entries(strchive: dict) -> list[dict]:
    """Flatten the merged catalog into the hg38 regions to fetch, one per locus."""
    entries = []
    for locus_id, meta in strchive.items():
        start, end = meta.get('start_hg38'), meta.get('stop_hg38')
        if not meta.get('chrom') or start is None or end is None:
            logger.warning(f'Skipping {locus_id}: no hg38 coordinates in the catalog')
            continue
        entries.append(
            {
                'chrom': meta['chrom'],
                'start': int(start),
                'end': int(end),
                'motifs': meta.get('reference_motif_reference_orientation', []),
                'locus_id': locus_id,
            }
        )
    return sorted(entries, key=lambda e: (e['chrom'], e['start']))


def _matches_locus(entry: dict, vcf_start: int, vcf_end: int) -> bool:
    """Check a variant's span against a BED entry within MATCH_TOLERANCE bp at both ends."""
    return abs(vcf_start - entry['start']) <= MATCH_TOLERANCE and abs(vcf_end - entry['end']) <= MATCH_TOLERANCE


def _format_value(variant, key: str, default: str = '.') -> str:
    """Read a single-sample FORMAT field as a display string, or default when absent."""
    try:
        values = variant.format(key)
    except KeyError:
        return default
    if values is None or len(values) == 0:
        return default
    first = values[0]
    # numeric fields come back as a per-sample array, strings as a flat array
    value = first[0] if hasattr(first, '__len__') and not isinstance(first, str) else first
    if isinstance(value, str):
        return value or default
    return f'{value:g}'


IUPAC_MAP = {
    'N': '[ACGT]',
    'R': '[AG]',
    'Y': '[CT]',
    'S': '[GC]',
    'W': '[AT]',
    'K': '[GT]',
    'M': '[AC]',
    'B': '[CGT]',
    'D': '[AGT]',
    'H': '[ACT]',
    'V': '[ACG]',
}

_motif_regex_cache: dict[str, re.Pattern] = {}


def _motif_to_regex(motif: str) -> re.Pattern:
    """Compile a motif with IUPAC ambiguity codes into a cached regex."""
    if motif not in _motif_regex_cache:
        pattern = ''.join(IUPAC_MAP.get(c, c) for c in motif.upper())
        _motif_regex_cache[motif] = re.compile(pattern)
    return _motif_regex_cache[motif]


def _is_degenerate(motif: str) -> bool:
    """Check whether a motif contains IUPAC ambiguity codes."""
    return any(c in IUPAC_MAP for c in motif.upper())


def count_motif_in_sequence(sequence: str, motif: str) -> int:
    """Count non-overlapping occurrences of a motif in a DNA sequence."""
    seq = sequence.upper()
    motif_upper = motif.upper()
    if not _is_degenerate(motif_upper):
        return seq.count(motif_upper)
    return len(_motif_to_regex(motif_upper).findall(seq))


def highlight_motifs_in_sequence(sequence: str, motif: str) -> str:
    """Wrap motif matches in HTML spans for colored highlighting."""
    seq = sequence.upper()
    motif_upper = motif.upper()
    pattern = _motif_to_regex(motif_upper)
    parts = []
    last_end = 0
    for m in pattern.finditer(seq):
        if m.start() > last_end:
            gap = seq[last_end : m.start()]
            parts.append(f'<span class="motif-interrupt">{gap}</span>')
        parts.append(f'<span class="motif-match">{m.group()}</span>')
        last_end = m.end()
    if last_end < len(seq):
        tail = seq[last_end:]
        parts.append(f'<span class="motif-interrupt">{tail}</span>')
    return ''.join(parts)


def compute_allele_repeat_units(
    gt_indices: list[int],
    alt_alleles: list[str],
    ref_seq: str,
    motif: str,
) -> tuple[float, float]:
    """Compute repeat unit counts for both alleles via motif counting."""
    ref_motif_count = count_motif_in_sequence(ref_seq, motif)

    alleles: list[float] = []
    for allele_idx in gt_indices:
        if 0 < allele_idx <= len(alt_alleles):
            alleles.append(float(count_motif_in_sequence(alt_alleles[allele_idx - 1], motif)))
        else:
            alleles.append(float(ref_motif_count))

    return alleles[0], alleles[1]


def classify_allele(repeat_units: float, locus_meta: dict) -> str:  # noqa: PLR0911
    """Classify a repeat count against STRchive thresholds, handling both expansion and contraction disorders."""
    benign_min = locus_meta.get('benign_min')
    benign_max = locus_meta.get('benign_max')
    intermediate_min = locus_meta.get('intermediate_min')
    intermediate_max = locus_meta.get('intermediate_max')
    pathogenic_min = locus_meta.get('pathogenic_min')
    pathogenic_max = locus_meta.get('pathogenic_max')

    # Contraction/discrete disorder: pathogenic range is below benign range.
    # Check benign first, then pathogenic, since ranges may overlap.
    if pathogenic_min is not None and benign_min is not None and pathogenic_min < benign_min:
        if benign_min <= repeat_units <= (benign_max if benign_max is not None else benign_min):
            return 'normal'
        if (
            intermediate_min is not None
            and intermediate_max is not None
            and intermediate_min <= repeat_units <= intermediate_max
        ):
            return 'intermediate'
        p_max = pathogenic_max if pathogenic_max is not None else pathogenic_min
        if pathogenic_min <= repeat_units <= p_max:
            return 'pathogenic'
        return 'uncertain'

    # Standard expansion logic
    if pathogenic_min is not None and repeat_units >= pathogenic_min:
        return 'pathogenic'
    if (
        intermediate_min is not None
        and intermediate_max is not None
        and intermediate_min <= repeat_units <= intermediate_max
    ):
        return 'intermediate'
    if benign_max is not None and repeat_units <= benign_max:
        return 'normal'
    if benign_max is not None and repeat_units > benign_max:
        if pathogenic_min is not None and repeat_units < pathogenic_min:
            return 'intermediate'
        return 'uncertain'
    return 'normal'


def classify_locus(status1: str, status2: str) -> str:
    """Return the more severe of two allele classifications."""
    return status1 if STATUS_PRIORITY.get(status1, 99) <= STATUS_PRIORITY.get(status2, 99) else status2


def _resolve_allele_seq(gt_idx: int, alt_alleles: list[str], ref_seq: str) -> str:
    """Return the allele sequence for a GT index (ref or alt)."""
    if 0 < gt_idx <= len(alt_alleles):
        return alt_alleles[gt_idx - 1]
    return ref_seq


def _parse_allreads(allreads_str: str, vcf_start: int, vcf_end: int, period: int) -> list[float]:
    """Parse the ALLREADS field into per-read repeat unit counts."""
    if not allreads_str or allreads_str == '.':
        return []
    ref_repeat_bp = vcf_end - vcf_start + 1
    ref_copies = ref_repeat_bp / period
    read_alleles: list[float] = []
    for token in allreads_str.split(';'):
        parts = token.strip().split('|')
        if not parts[0]:
            continue
        try:
            bp_diff = int(parts[0])
            count = int(parts[1]) if len(parts) > 1 else 1
            ru = math.floor(ref_copies + bp_diff / period)
            read_alleles.extend([ru] * count)
        except ValueError:
            pass
    return read_alleles


def _build_external_links(meta: dict, locus_id: str) -> list[dict]:
    """Build external database links from STRchive metadata."""
    links = [{'label': 'STRchive', 'url': f'https://strchive.org/loci/{locus_id.lower()}/'}]
    for key, label, url_template in EXTERNAL_LINK_DEFS:
        ids = meta.get(key, [])
        if ids:
            links.append({'label': label, 'url': url_template.format(ids[0])})
    return links


def _threshold_parts(meta: dict) -> list[str]:
    """Summarise a locus's thresholds for display, enumerating counts for discrete/contraction loci."""
    benign_min = meta.get('benign_min')
    benign_max = meta.get('benign_max')
    intermediate_min = meta.get('intermediate_min')
    intermediate_max = meta.get('intermediate_max')
    pathogenic_min = meta.get('pathogenic_min')
    pathogenic_max = meta.get('pathogenic_max')

    # Discrete/contraction loci (VWA1, MIR7-2) put their pathogenic range below the benign one,
    # so '<=benign / >=pathogenic' would contradict itself. The ranges are only a few counts
    # wide, so label each one via classify_allele rather than restating its precedence rules.
    if pathogenic_min is not None and benign_min is not None and pathogenic_min < benign_min:
        lo = int(pathogenic_min)
        hi = int(max(benign_max or benign_min, pathogenic_max or pathogenic_min))
        by_status: dict[str, list[int]] = defaultdict(list)
        for count in range(lo, hi + 1):
            by_status[classify_allele(count, meta)].append(count)

        parts = [
            f'{label}: ' + ', '.join(str(c) for c in by_status[status])
            for status, label in (('normal', 'Normal'), ('intermediate', 'Intermediate'), ('pathogenic', 'Pathogenic'))
            if by_status[status]
        ]
        parts.append('any other count: uncertain')
        return parts

    parts = []
    if benign_max is not None:
        parts.append(f'Normal: ≤{benign_max}')
    if intermediate_min is not None and intermediate_max is not None:
        parts.append(f'Intermediate: {intermediate_min}-{intermediate_max}')
    if pathogenic_min is not None:
        parts.append(f'Pathogenic: ≥{pathogenic_min}')
    return parts


def _build_locus_meta(meta: dict, entry: dict) -> dict:
    """Assemble shared locus metadata from STRchive and BED entry fields."""
    motif_list = meta.get('reference_motif_reference_orientation', entry['motifs'])
    return {
        'gene': meta.get('gene', ''),
        'disease': meta.get('disease', ''),
        'disease_id': meta.get('disease_id', ''),
        'chrom': entry['chrom'],
        'start': entry['start'],
        'end': entry['end'],
        'motif': ','.join(motif_list),
        'primary_motif': motif_list[0] if motif_list else '',
        'period': meta.get('motif_len', 0),
        'location_in_gene': meta.get('location_in_gene', ''),
        'inheritance': ', '.join(meta.get('inheritance', [])),
        'mechanism': meta.get('mechanism', ''),
        'benign_min': meta.get('benign_min'),
        'benign_max': meta.get('benign_max'),
        'intermediate_min': meta.get('intermediate_min'),
        'intermediate_max': meta.get('intermediate_max'),
        'pathogenic_min': meta.get('pathogenic_min'),
        'pathogenic_max': meta.get('pathogenic_max'),
        'ref_copies': meta.get('ref_copies'),
        'thresholds': _threshold_parts(meta),
        'evidence': ', '.join(meta.get('evidence', [])),
        'external_links': _build_external_links(meta, entry['locus_id']),
    }


def scan_vcf(vcf_path: str, entries: list[dict], strchive: dict, sex: str = 'unknown') -> list[dict]:
    """Fetch each disease locus from an indexed VCF, returning results sorted worst-first."""
    results = {}
    vcf = VCF(vcf_path)

    for entry in entries:
        locus_id = entry['locus_id']
        meta = strchive[locus_id]

        # Pad the fetch by the tolerance, since LongTR's coordinates do not exactly match
        # STRchive's; the same fuzzy comparison then filters whatever the index returns.
        start = max(1, entry['start'] - MATCH_TOLERANCE)
        region = f'{entry["chrom"]}:{start}-{entry["end"] + MATCH_TOLERANCE}'
        for variant in vcf(region):
            vcf_start = variant.INFO.get('START') or variant.POS
            vcf_end = variant.INFO.get('END') or vcf_start
            if not _matches_locus(entry, vcf_start, vcf_end):
                continue
            results[locus_id] = _process_vcf_record(variant, entry, meta, vcf_start, vcf_end, sex=sex)
            break

    vcf.close()
    _add_missing_loci(results, entries, strchive)

    return sorted(
        results.values(),
        key=lambda r: (
            STATUS_PRIORITY.get(r['locus_status'], 99),
            r['chrom'],
            r['start'],
        ),
    )


def _process_vcf_record(variant, match, meta, vcf_start, vcf_end, *, sex: str = 'unknown') -> dict:
    """Extract genotype data and classify alleles for a single VCF record."""
    period = variant.INFO.get('PERIOD') or meta.get('motif_len', 3)
    alt_alleles = list(variant.ALT)

    # cyvcf2 gives [allele1, ..., alleleN, phased] - the phase flag is always last, so slice it
    # off rather than taking a fixed two. A haploid call yields one index, duplicated below so
    # downstream code can keep assuming a pair; -1 marks a missing call, which counts as ref.
    gt = variant.genotypes[0] if variant.genotypes else [0, 0, False]
    gt_indices = [max(0, i) for i in gt[:-1]]
    if len(gt_indices) == 1:
        gt_indices *= 2

    # Shared locus metadata from STRchive + BED, with VCF-specific overrides
    base = _build_locus_meta(meta, match)
    base['start'] = vcf_start
    base['end'] = vcf_end
    base['period'] = period

    ref_seq = variant.REF
    first, second = compute_allele_repeat_units(
        gt_indices,
        alt_alleles,
        ref_seq,
        base['primary_motif'],
    )
    seq_first = _resolve_allele_seq(gt_indices[0], alt_alleles, ref_seq)
    seq_second = _resolve_allele_seq(gt_indices[1], alt_alleles, ref_seq)

    # Males carry one X, and every chrX disease locus in STRchive sits outside the
    # pseudoautosomal regions. LongTR calls chrX diploid regardless of sex, so collapse to a
    # single allele: all chrX loci are expansion disorders, making the larger count the
    # conservative choice (never under-call an expansion).
    hemizygous = sex == 'male' and match['chrom'] == 'chrX'
    if hemizygous and second > first:
        first, seq_first = second, seq_second

    a1 = first
    seq1 = seq_first
    a2: float | None = None if hemizygous else second
    seq2: str | None = None if hemizygous else seq_second

    s1 = classify_allele(a1, meta)
    s2: str | None = None if a2 is None else classify_allele(a2, meta)

    locus_status = s1 if s2 is None else classify_locus(s1, s2)

    # A female heterozygous for a pathogenic allele at an X-linked *recessive* locus is not
    # straightforwardly affected - her second X compensates - so the call is downgraded to
    # uncertain rather than asserting disease. Males are hemizygous and keep the pathogenic
    # call; X-linked dominant loci (FMR1) manifest in females and are left alone.
    if (
        sex == 'female'
        and match['chrom'] == 'chrX'
        and 'XR' in meta.get('inheritance', [])
        and [s1, s2].count('pathogenic') == 1
    ):
        locus_status = 'uncertain'

    base.update(
        {
            'locus_id': match['locus_id'],
            'allele1_ru': a1,
            'allele2_ru': a2,
            'allele1_seq': seq1,
            'allele2_seq': seq2,
            'allele1_status': s1,
            'allele2_status': s2,
            'hemizygous': hemizygous,
            'locus_status': locus_status,
            'dp': _format_value(variant, 'DP'),
            'q': _format_value(variant, 'Q'),
            'pq': _format_value(variant, 'PQ'),
            'gldiff': _format_value(variant, 'GLDIFF'),
            'gt': ('|' if gt[-1] else '/').join('.' if i < 0 else str(i) for i in gt[:-1]),
            'filter': variant.FILTER or 'PASS',
            'read_alleles': _parse_allreads(_format_value(variant, 'ALLREADS', ''), vcf_start, vcf_end, period),
            'genotyped': True,
        }
    )
    return base


def _add_missing_loci(results: dict, entries: list[dict], strchive: dict) -> None:
    """Add not-genotyped placeholder entries for loci absent from the VCF."""
    for entry in entries:
        lid = entry['locus_id']
        if lid in results:
            continue
        meta = strchive[lid]
        base = _build_locus_meta(meta, entry)
        base.update(
            {
                'locus_id': lid,
                'allele1_ru': None,
                'allele2_ru': None,
                'allele1_seq': None,
                'allele2_seq': None,
                'allele1_status': 'not_genotyped',
                'allele2_status': 'not_genotyped',
                'hemizygous': False,
                'locus_status': 'not_genotyped',
                'dp': '.',
                'q': '.',
                'pq': '.',
                'gldiff': '.',
                'gt': '.',
                'filter': '.',
                'read_alleles': [],
                'genotyped': False,
            }
        )
        results[lid] = base


def _compute_gauge_model(result: dict, width: int = 600) -> dict | None:
    """Compute gauge geometry: zones, ticks, and allele markers as coordinates."""
    if not result['genotyped']:
        return None

    # Discrete/contraction loci span only a few counts with the pathogenic range below the
    # benign one, so a continuous bar would paint overlapping zones that contradict the
    # enumerated threshold text. Those few values are better read from the text alone.
    benign_min = result.get('benign_min')
    pathogenic_min = result.get('pathogenic_min')
    if pathogenic_min is not None and benign_min is not None and pathogenic_min < benign_min:
        return None

    a1 = result['allele1_ru']
    a2 = result['allele2_ru']
    scale_values = [
        v for v in [a1, a2, result['benign_max'], result['pathogenic_min'], result['intermediate_max']] if v is not None
    ]
    if not scale_values:
        return None

    scale_max = max(max(scale_values) * 1.3, 10)
    h, bar_y, bar_h, margin_l = 70, 25, 20, 10
    bar_w = width - 2 * margin_l

    def x_pos(val) -> float:
        return margin_l + (val / scale_max) * bar_w

    benign_max = result['benign_max']
    pathogenic_min = result['pathogenic_min']
    intermediate_min = result['intermediate_min']
    intermediate_max = result['intermediate_max']

    zones = []
    if benign_max is not None:
        zones.append({'x': margin_l, 'w': x_pos(benign_max) - margin_l, 'color': '#28a745', 'rx': 3})
    if intermediate_min is not None and intermediate_max is not None:
        ix1, ix2 = x_pos(intermediate_min), x_pos(intermediate_max)
        zones.append({'x': ix1, 'w': ix2 - ix1, 'color': '#ffc107'})
    elif benign_max is not None and pathogenic_min is not None and pathogenic_min > benign_max + 1:
        ix1, ix2 = x_pos(benign_max), x_pos(pathogenic_min)
        zones.append({'x': ix1, 'w': ix2 - ix1, 'color': '#ffc107'})
    if pathogenic_min is not None:
        px1 = x_pos(pathogenic_min)
        zones.append({'x': px1, 'w': margin_l + bar_w - px1, 'color': '#dc3545', 'rx': 3})

    tick_values = sorted(
        {v for v in [0, benign_max, intermediate_min, intermediate_max, pathogenic_min] if v is not None},
    )
    min_gap = 30
    last_x = -999.0
    ticks = []
    for tv in tick_values:
        tx = x_pos(tv)
        show_label = tx - last_x >= min_gap
        ticks.append({'x': tx, 'label': f'{tv:.0f}', 'show_label': show_label})
        if show_label:
            last_x = tx

    read_dots = [{'cx': x_pos(ra)} for ra in result.get('read_alleles', []) if ra <= scale_max]

    markers = []
    for allele_val, color in [(a1, '#0d6efd'), (a2, '#6610f2')]:
        if allele_val is None:
            continue
        ax = x_pos(min(allele_val, scale_max))
        label = f'{allele_val:.0f}' if allele_val == int(allele_val) else f'{allele_val:.1f}'
        markers.append({'x': ax, 'color': color, 'label': label})

    return {
        'width': width,
        'h': h,
        'bar_y': bar_y,
        'bar_h': bar_h,
        'margin_l': margin_l,
        'bar_w': bar_w,
        'zones': zones,
        'ticks': ticks,
        'read_dots': read_dots,
        'markers': markers,
    }


def status_badge(status: str) -> str:
    """Return an HTML badge span for a classification status."""
    colors = {
        'pathogenic': ('#dc3545', '#fff'),
        'intermediate': ('#ffc107', '#333'),
        'uncertain': ('#fd7e14', '#fff'),
        'normal': ('#28a745', '#fff'),
        'not_genotyped': ('#6c757d', '#fff'),
    }
    bg, fg = colors.get(status, ('#6c757d', '#fff'))
    label = status.replace('_', ' ').title()
    return f'<span class="badge" style="background:{bg};color:{fg}">{label}</span>'


def _summarise_results(results: list[dict]) -> dict[str, int]:
    """Count results by genotyped status and locus classification."""
    counts: dict[str, int] = {'total_loci': len(results), 'genotyped': 0, 'not_genotyped': 0}
    for r in results:
        if r['genotyped']:
            counts['genotyped'] += 1
            counts[r['locus_status']] = counts.get(r['locus_status'], 0) + 1
        else:
            counts['not_genotyped'] += 1
    return counts


def _build_participant(
    sample_id: str,
    sex: str,
    birth_year: str,
    age_of_onset: str,
    hpo_terms: str,
) -> dict:
    """Assemble the participant header fields, deriving age from the recorded birth year."""
    age = None
    if birth_year.isdigit():
        age = datetime.datetime.now(tz=datetime.timezone.utc).year - int(birth_year)
    return {
        'sample_id': sample_id,
        'sex': sex,
        'age': age,
        'age_of_onset': age_of_onset,
        'hpo_terms': [t.strip() for t in hpo_terms.split(',') if t.strip()],
    }


def generate_html(
    results: list[dict],
    sample_name: str,
    summary: dict[str, int],
    report_type: str,
    participant: dict | None = None,
) -> str:
    """Render the full HTML report from results via the Jinja2 template."""
    template_dir = Path(__file__).resolve().parent / 'templates'
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(template_dir)),
        autoescape=True,
    )
    env.globals['gauge_model'] = _compute_gauge_model
    env.globals['status_badge'] = lambda s: Markup(status_badge(s))  # noqa: S704
    env.globals['highlight_seq'] = lambda seq, motif: Markup(highlight_motifs_in_sequence(seq, motif))  # noqa: S704

    template = env.get_template('longtr_pathogenic.html.jinja')

    evidence_levels = sorted({r.get('evidence', '') for r in results if r.get('evidence')})
    evidence_counts = {level: sum(1 for r in results if r.get('evidence') == level) for level in evidence_levels}

    # A locus can carry several modes (e.g. 'AD, AR'), so count membership rather than the whole string
    modes_per_result = [[m.strip() for m in r.get('inheritance', '').split(',') if m.strip()] for r in results]
    inheritance_levels = sorted({m for modes in modes_per_result for m in modes})
    inheritance_counts = {mode: sum(1 for modes in modes_per_result if mode in modes) for mode in inheritance_levels}

    return template.render(
        sample_name=sample_name,
        participant=participant or {},
        report_type=report_type,
        results=results,
        n_genotyped=summary.get('genotyped', 0),
        n_pathogenic=summary.get('pathogenic', 0),
        n_intermediate=summary.get('intermediate', 0),
        n_uncertain=summary.get('uncertain', 0),
        n_normal=summary.get('normal', 0),
        n_missing=summary.get('not_genotyped', 0),
        evidence_levels=evidence_levels,
        evidence_counts=evidence_counts,
        inheritance_levels=inheritance_levels,
        inheritance_counts=inheritance_counts,
        inheritance_names=INHERITANCE_NAMES,
    )


def build_json_output(
    results: list[dict],
    sample_name: str,
    summary: dict[str, int],
    strchive_json_path: str,
    custom_loci_json_path: str,
) -> dict:
    """Build the structured JSON output dict from screening results."""
    json_fields = (
        'locus_id',
        'gene',
        'disease',
        'disease_id',
        'chrom',
        'start',
        'end',
        'motif',
        'primary_motif',
        'period',
        'inheritance',
        'mechanism',
        'allele1_ru',
        'allele2_ru',
        'allele1_seq',
        'allele2_seq',
        'allele1_status',
        'allele2_status',
        'hemizygous',
        'locus_status',
        'benign_min',
        'benign_max',
        'intermediate_min',
        'intermediate_max',
        'pathogenic_min',
        'pathogenic_max',
        'ref_copies',
        'thresholds',
        'gt',
        'dp',
        'q',
        'evidence',
        'external_links',
        'genotyped',
    )
    return {
        'sample_name': sample_name,
        'catalog': {
            'strchive_json': strchive_json_path,
            'custom_loci_json': custom_loci_json_path,
        },
        'summary': summary,
        'loci': [{k: v for k, v in r.items() if k in json_fields} for r in results],
    }


def generate_report(
    vcf_path: str,
    strchive_json: str,
    output_html: str,
    output_json: str,
    *,
    report_type: str,
    sample_id: str,
    loci_list: set[str] | None = None,
    sex: str = 'unknown',
    birth_year: str = '',
    age_of_onset: str = '',
    hpo_terms: str = '',
    custom_loci_json: str = '',
):
    """Load references, scan VCF, optionally filter by loci list, and write outputs."""
    strchive = load_strchive_json(strchive_json, custom_loci_json or None)
    results = scan_vcf(vcf_path, catalog_entries(strchive), strchive, sex)

    if loci_list:
        results = [r for r in results if r['locus_id'] in loci_list]

    summary = _summarise_results(results)

    participant = _build_participant(sample_id, sex, birth_year, age_of_onset, hpo_terms)
    html_content = generate_html(results, sample_id, summary, report_type, participant)
    with open(output_html, 'w') as f:
        f.write(html_content)

    json_output = build_json_output(results, sample_id, summary, strchive_json, custom_loci_json)
    with open(output_json, 'w') as f:
        json.dump(json_output, f, indent=2)

    def _fmt_ru(v) -> str:
        return f'{v:.0f}' if v == int(v) else f'{v:.1f}'

    logger.info(f'Screened {len(results)} disease loci ({summary["genotyped"]} genotyped)')
    for r in results:
        if r['locus_status'] in ('pathogenic', 'intermediate', 'uncertain'):
            counts = _fmt_ru(r['allele1_ru'])
            if r['allele2_ru'] is not None:
                counts += f'/{_fmt_ru(r["allele2_ru"])}'
            elif r.get('hemizygous'):
                counts += ' (hemizygous)'
            logger.warning(f'{r["gene"]} ({r["disease"]}): {r["locus_status"]} — {counts} repeats')
    logger.info(f'HTML report: {output_html}')
    logger.info(f'JSON results: {output_json}')


def cli_main():
    """Parse CLI arguments and run the report pipeline."""
    parser = ArgumentParser(
        description='Screen a LongTR VCF against STRchive disease-associated TR loci.',
    )
    parser.add_argument('--vcf_path', required=True, help='Path to LongTR VCF file')
    parser.add_argument('--strchive_json', required=True, help='Path to STRchive-loci.json (upstream)')
    parser.add_argument(
        '--custom_loci_json',
        default='',
        help='Path to locally curated loci overlaid on the STRchive catalog',
    )
    parser.add_argument('--output_html', required=True, help='Output HTML file')
    parser.add_argument('--output_json', required=True, help='Output JSON file')
    parser.add_argument('--report_type', required=True, help='Loci list name this report covers')
    parser.add_argument('--loci_list', help='Locus IDs to include', nargs='+')
    parser.add_argument(
        '--sex',
        default='unknown',
        choices=['male', 'female', 'unknown'],
        help='Reported sex; males are treated as hemizygous at chrX loci',
    )
    parser.add_argument('--sample_id', required=True, help='Sample ID to display on the report')
    parser.add_argument('--birth_year', default='', help='Participant birth year, if recorded in metamist')
    parser.add_argument('--age_of_onset', default='', help='Participant age of onset, if recorded in metamist')
    parser.add_argument('--hpo_terms', default='', help='Comma-separated HPO term IDs, if recorded in metamist')
    args = parser.parse_args()

    loci_set = set(args.loci_list) if args.loci_list else None

    generate_report(
        args.vcf_path,
        args.strchive_json,
        args.output_html,
        args.output_json,
        report_type=args.report_type,
        sample_id=args.sample_id,
        loci_list=loci_set,
        sex=args.sex,
        birth_year=args.birth_year,
        age_of_onset=args.age_of_onset,
        hpo_terms=args.hpo_terms,
        custom_loci_json=args.custom_loci_json,
    )


if __name__ == '__main__':
    cli_main()
