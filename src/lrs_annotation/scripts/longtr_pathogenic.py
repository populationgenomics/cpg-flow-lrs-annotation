"""
Screen a LongTR VCF against STRchive disease-associated TR loci and
generate a standalone HTML report with STRIPY-style gauge visualizations.

Requires two reference files from STRchive (github.com/dashnowlab/STRchive):
  - STRchive-loci.json (disease locus metadata + thresholds)
  - STRchive-disease-loci.hg38.longTR.bed (coordinates for matching)
"""

import gzip
import json
import re
import urllib.request
from argparse import ArgumentParser
from collections import defaultdict
from importlib import resources
from pathlib import Path

import jinja2
from markupsafe import Markup

STRCHIVE_JSON_URL = (
    'https://raw.githubusercontent.com/dashnowlab/STRchive/main/data/STRchive-loci.json'
)
STRCHIVE_BED_URL = (
    'https://raw.githubusercontent.com/dashnowlab/STRchive/main/data/catalogs/'
    'STRchive-disease-loci.hg38.longTR.bed'
)

MATCH_TOLERANCE = 20


def open_vcf(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path)


def download_if_missing(path: Path, url: str) -> Path:
    if path.exists():
        return path
    print(f'Downloading {url} ...')
    urllib.request.urlretrieve(url, path)
    return path


def load_strchive_json(path: str) -> dict:
    with open(path) as f:
        loci = json.load(f)
    return {locus['id']: locus for locus in loci}


def load_longtr_bed(path: str) -> list[dict]:
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            entries.append({
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'motifs': parts[3].split(','),
                'locus_id': parts[4],
            })
    return entries


def build_locus_index(bed_entries: list[dict]) -> dict:
    index = defaultdict(list)
    for entry in bed_entries:
        index[entry['chrom']].append(entry)
    return dict(index)


def find_matching_locus(chrom: str, vcf_start: int, vcf_end: int, index: dict) -> dict | None:
    candidates = index.get(chrom, [])
    for entry in candidates:
        if abs(vcf_start - entry['start']) <= MATCH_TOLERANCE and abs(vcf_end - entry['end']) <= MATCH_TOLERANCE:
            return entry
    return None


def parse_info(info_str: str) -> dict[str, str]:
    fields = {}
    for item in info_str.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            fields[k] = v
        else:
            fields[item] = 'true'
    return fields


def parse_format_sample(fmt_str: str, sample_str: str) -> dict[str, str]:
    keys = fmt_str.split(':')
    values = sample_str.split(':')
    return dict(zip(keys, values, strict=False))


IUPAC_MAP = {
    'N': '[ACGT]', 'R': '[AG]', 'Y': '[CT]', 'S': '[GC]',
    'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
    'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
}

_motif_regex_cache: dict[str, re.Pattern] = {}


def _motif_to_regex(motif: str) -> re.Pattern:
    if motif not in _motif_regex_cache:
        pattern = ''.join(IUPAC_MAP.get(c, c) for c in motif.upper())
        _motif_regex_cache[motif] = re.compile(pattern)
    return _motif_regex_cache[motif]


def _is_degenerate(motif: str) -> bool:
    return any(c in IUPAC_MAP for c in motif.upper())


def count_motif_in_sequence(sequence: str, motif: str) -> int:
    seq = sequence.upper()
    motif_upper = motif.upper()
    if not _is_degenerate(motif_upper):
        return seq.count(motif_upper)
    return len(_motif_to_regex(motif_upper).findall(seq))


def highlight_motifs_in_sequence(sequence: str, motif: str) -> str:
    seq = sequence.upper()
    motif_upper = motif.upper()
    pattern = _motif_to_regex(motif_upper)
    parts = []
    last_end = 0
    for m in pattern.finditer(seq):
        if m.start() > last_end:
            gap = seq[last_end:m.start()]
            parts.append(f'<span class="motif-interrupt">{gap}</span>')
        parts.append(f'<span class="motif-match">{m.group()}</span>')
        last_end = m.end()
    if last_end < len(seq):
        tail = seq[last_end:]
        parts.append(f'<span class="motif-interrupt">{tail}</span>')
    return ''.join(parts)


def compute_allele_repeat_units(
    ref_seq: str,
    period: int,
    gb_str: str,
    info_start: int,
    info_end: int,
    gt_str: str,
    alt_alleles: list[str],
    motif: str,
) -> tuple[float, float]:
    ref_repeat_bp = info_end - info_start + 1
    ref_copies = ref_repeat_bp / period

    sep = '|' if '|' in gt_str else '/'
    gt_indices = []
    for idx_str in gt_str.split(sep):
        try:
            gt_indices.append(int(idx_str))
        except ValueError:
            gt_indices.append(0)
    if len(gt_indices) == 1:
        gt_indices = [gt_indices[0], gt_indices[0]]

    gb_sep = '|' if '|' in gb_str else '/'
    gb_parts = gb_str.split(gb_sep)

    alleles = []
    for i, allele_idx in enumerate(gt_indices):
        if allele_idx > 0 and allele_idx <= len(alt_alleles) and alt_alleles[allele_idx - 1] != '.':
            alleles.append(count_motif_in_sequence(alt_alleles[allele_idx - 1], motif))
        else:
            try:
                bp_diff = int(gb_parts[i]) if i < len(gb_parts) else 0
            except ValueError:
                bp_diff = 0
            alleles.append(round(ref_copies + bp_diff / period, 1))

    return alleles[0], alleles[1]


def classify_allele(repeat_units: float, locus_meta: dict) -> str:
    benign_max = locus_meta.get('benign_max')
    intermediate_min = locus_meta.get('intermediate_min')
    intermediate_max = locus_meta.get('intermediate_max')
    pathogenic_min = locus_meta.get('pathogenic_min')

    if pathogenic_min is not None and repeat_units >= pathogenic_min:
        return 'pathogenic'
    if intermediate_min is not None and intermediate_max is not None:
        if intermediate_min <= repeat_units <= intermediate_max:
            return 'intermediate'
    if benign_max is not None and repeat_units <= benign_max:
        return 'normal'
    if benign_max is not None and pathogenic_min is not None:
        if repeat_units > benign_max and repeat_units < pathogenic_min:
            return 'intermediate'
    if benign_max is not None and repeat_units > benign_max:
        return 'uncertain'
    return 'normal'


def classify_locus(status1: str, status2: str) -> str:
    priority = {'pathogenic': 0, 'intermediate': 1, 'uncertain': 2, 'normal': 3}
    return status1 if priority.get(status1, 99) <= priority.get(status2, 99) else status2


def scan_vcf(vcf_path: str, locus_index: dict, strchive: dict) -> list[dict]:
    results = {}
    sample_name = ''

    with open_vcf(vcf_path) as f:
        for line in f:
            if isinstance(line, bytes):
                line = line.decode()
            line = line.rstrip('\n')

            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                cols = line.split('\t')
                if len(cols) >= 10:
                    sample_name = cols[9]
                continue

            cols = line.split('\t')
            if len(cols) < 10:
                continue

            chrom = cols[0]
            if chrom not in locus_index:
                continue

            info = parse_info(cols[7])
            vcf_start = int(info.get('START', cols[1]))
            vcf_end = int(info.get('END', vcf_start))

            match = find_matching_locus(chrom, vcf_start, vcf_end, locus_index)
            if match is None:
                continue

            locus_id = match['locus_id']
            if locus_id in results:
                continue

            meta = strchive.get(locus_id)
            if meta is None:
                continue

            sample_data = parse_format_sample(cols[8], cols[9])
            period = int(info.get('PERIOD', meta.get('motif_len', 3)))

            alt_field = cols[4]
            alt_alleles = alt_field.split(',') if alt_field != '.' else []

            motif_list = meta.get('reference_motif_reference_orientation', match['motifs'])
            primary_motif = motif_list[0] if motif_list else match['motifs'][0]

            gt_str = sample_data.get('GT', '0/0')
            gb = sample_data.get('GB', '0|0')
            a1, a2 = compute_allele_repeat_units(
                cols[3], period, gb, vcf_start, vcf_end,
                gt_str, alt_alleles, primary_motif,
            )

            gt_sep = '|' if '|' in gt_str else '/'
            gt_indices = []
            for idx_str in gt_str.split(gt_sep):
                try:
                    gt_indices.append(int(idx_str))
                except ValueError:
                    gt_indices.append(0)
            if len(gt_indices) == 1:
                gt_indices = [gt_indices[0], gt_indices[0]]

            ref_seq = cols[3]
            allele1_seq = alt_alleles[gt_indices[0] - 1] if gt_indices[0] > 0 and gt_indices[0] <= len(alt_alleles) else ref_seq
            allele2_seq = alt_alleles[gt_indices[1] - 1] if gt_indices[1] > 0 and gt_indices[1] <= len(alt_alleles) else ref_seq

            s1 = classify_allele(a1, meta)
            s2 = classify_allele(a2, meta)

            allreads_str = sample_data.get('ALLREADS', '')
            read_alleles = []
            if allreads_str and allreads_str != '.':
                ref_repeat_bp = vcf_end - vcf_start + 1
                ref_copies = ref_repeat_bp / period
                for entry in allreads_str.split(';'):
                    entry = entry.strip()
                    if not entry:
                        continue
                    parts = entry.split('|')
                    if len(parts) >= 2:
                        try:
                            bp_diff = int(parts[0])
                            count = int(parts[1])
                            ru = round(ref_copies + bp_diff / period, 1)
                            for _ in range(count):
                                read_alleles.append(ru)
                        except ValueError:
                            pass

            results[locus_id] = {
                'locus_id': locus_id,
                'gene': meta.get('gene', ''),
                'disease': meta.get('disease', ''),
                'disease_id': meta.get('disease_id', ''),
                'chrom': chrom,
                'start': vcf_start,
                'end': vcf_end,
                'motif': ','.join(motif_list),
                'primary_motif': primary_motif,
                'period': period,
                'location_in_gene': meta.get('location_in_gene', ''),
                'inheritance': ', '.join(meta.get('inheritance', [])),
                'mechanism': meta.get('mechanism', ''),
                'allele1_ru': a1,
                'allele2_ru': a2,
                'allele1_seq': allele1_seq,
                'allele2_seq': allele2_seq,
                'allele1_status': s1,
                'allele2_status': s2,
                'locus_status': classify_locus(s1, s2),
                'benign_min': meta.get('benign_min'),
                'benign_max': meta.get('benign_max'),
                'intermediate_min': meta.get('intermediate_min'),
                'intermediate_max': meta.get('intermediate_max'),
                'pathogenic_min': meta.get('pathogenic_min'),
                'pathogenic_max': meta.get('pathogenic_max'),
                'ref_copies': meta.get('ref_copies'),
                'dp': sample_data.get('DP', '.'),
                'q': sample_data.get('Q', '.'),
                'pq': sample_data.get('PQ', '.'),
                'gldiff': sample_data.get('GLDIFF', '.'),
                'gt': sample_data.get('GT', '.'),
                'filter': sample_data.get('FILTER', '.'),
                'read_alleles': read_alleles,
                'genotyped': True,
            }

    for entry in sum(locus_index.values(), []):
        lid = entry['locus_id']
        if lid not in results:
            meta = strchive.get(lid, {})
            motif_list = meta.get('reference_motif_reference_orientation', entry['motifs'])
            results[lid] = {
                'locus_id': lid,
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
                'allele1_ru': None,
                'allele2_ru': None,
                'allele1_seq': None,
                'allele2_seq': None,
                'allele1_status': 'not_genotyped',
                'allele2_status': 'not_genotyped',
                'locus_status': 'not_genotyped',
                'benign_min': meta.get('benign_min'),
                'benign_max': meta.get('benign_max'),
                'intermediate_min': meta.get('intermediate_min'),
                'intermediate_max': meta.get('intermediate_max'),
                'pathogenic_min': meta.get('pathogenic_min'),
                'pathogenic_max': meta.get('pathogenic_max'),
                'ref_copies': meta.get('ref_copies'),
                'dp': '.',
                'q': '.',
                'pq': '.',
                'gldiff': '.',
                'gt': '.',
                'filter': '.',
                'read_alleles': [],
                'genotyped': False,
            }

    sorted_results = sorted(results.values(), key=lambda r: (
        {'pathogenic': 0, 'intermediate': 1, 'uncertain': 2, 'normal': 3, 'not_genotyped': 4}.get(r['locus_status'], 5),
        r['chrom'],
        r['start'],
    ))

    return sorted_results, sample_name


def svg_gauge(result: dict, width: int = 600) -> str:
    if not result['genotyped']:
        return '<div class="gauge-placeholder">Not genotyped — locus not found in VCF</div>'

    benign_max = result['benign_max']
    pathogenic_min = result['pathogenic_min']
    pathogenic_max = result['pathogenic_max']
    intermediate_min = result['intermediate_min']
    intermediate_max = result['intermediate_max']
    a1 = result['allele1_ru']
    a2 = result['allele2_ru']

    scale_values = [v for v in [a1, a2, benign_max, pathogenic_min, intermediate_max] if v is not None]
    if not scale_values:
        return '<div class="gauge-placeholder">No threshold data available</div>'

    scale_max = max(scale_values) * 1.3
    scale_max = max(scale_max, 10)

    h = 70
    bar_y = 25
    bar_h = 20
    margin_l = 10
    bar_w = width - 2 * margin_l

    def x_pos(val):
        return margin_l + (val / scale_max) * bar_w

    svg_parts = [f'<svg viewBox="0 0 {width} {h}" xmlns="http://www.w3.org/2000/svg" class="gauge-svg">']

    svg_parts.append(f'<rect x="{margin_l}" y="{bar_y}" width="{bar_w}" height="{bar_h}" fill="#e9ecef" rx="3"/>')

    if benign_max is not None:
        bw = x_pos(min(benign_max, scale_max)) - margin_l
        svg_parts.append(f'<rect x="{margin_l}" y="{bar_y}" width="{bw}" height="{bar_h}" fill="#28a745" opacity="0.5" rx="3"/>')

    if intermediate_min is not None and intermediate_max is not None:
        ix1 = x_pos(intermediate_min)
        ix2 = x_pos(min(intermediate_max, scale_max))
        svg_parts.append(f'<rect x="{ix1}" y="{bar_y}" width="{ix2-ix1}" height="{bar_h}" fill="#ffc107" opacity="0.5"/>')
    elif benign_max is not None and pathogenic_min is not None and pathogenic_min > benign_max + 1:
        ix1 = x_pos(benign_max)
        ix2 = x_pos(pathogenic_min)
        svg_parts.append(f'<rect x="{ix1}" y="{bar_y}" width="{ix2-ix1}" height="{bar_h}" fill="#ffc107" opacity="0.5"/>')

    if pathogenic_min is not None:
        px1 = x_pos(pathogenic_min)
        px2 = margin_l + bar_w
        svg_parts.append(f'<rect x="{px1}" y="{bar_y}" width="{px2-px1}" height="{bar_h}" fill="#dc3545" opacity="0.5" rx="3"/>')

    tick_values = set()
    for v in [0, benign_max, intermediate_min, intermediate_max, pathogenic_min]:
        if v is not None:
            tick_values.add(v)
    min_label_gap = 30
    last_label_x = -999
    for tv in sorted(tick_values):
        tx = x_pos(tv)
        svg_parts.append(f'<line x1="{tx}" y1="{bar_y}" x2="{tx}" y2="{bar_y+bar_h}" stroke="#333" stroke-width="0.5" stroke-dasharray="2,1"/>')
        if tx - last_label_x >= min_label_gap:
            svg_parts.append(f'<text x="{tx}" y="{bar_y+bar_h+11}" text-anchor="middle" font-size="9" fill="#666">{tv:.0f}</text>')
            last_label_x = tx

    read_alleles = result.get('read_alleles', [])
    if read_alleles:
        for ra in read_alleles:
            if ra <= scale_max:
                rx = x_pos(ra)
                svg_parts.append(f'<circle cx="{rx}" cy="{bar_y + bar_h/2}" r="2" fill="rgba(100,100,100,0.15)" stroke="none"/>')

    for allele_val, label, color in [(a1, 'A1', '#0d6efd'), (a2, 'A2', '#6610f2')]:
        if allele_val is None:
            continue
        ax = x_pos(min(allele_val, scale_max))
        svg_parts.append(f'<line x1="{ax}" y1="{bar_y-2}" x2="{ax}" y2="{bar_y+bar_h+2}" stroke="{color}" stroke-width="2.5"/>')
        svg_parts.append(f'<circle cx="{ax}" cy="{bar_y-5}" r="4" fill="{color}"/>')
        label_text = f'{allele_val:.0f}' if allele_val == int(allele_val) else f'{allele_val:.1f}'
        svg_parts.append(f'<text x="{ax}" y="{bar_y-10}" text-anchor="middle" font-size="8" font-weight="bold" fill="{color}">{label_text}</text>')

    svg_parts.append('</svg>')
    return '\n'.join(svg_parts)


def status_badge(status: str) -> str:
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


def generate_html(results: list[dict], sample_name: str) -> str:
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(
            str(resources.files('lrs_annotation') / 'templates'),
        ),
        autoescape=True,
    )
    env.globals['svg_gauge'] = lambda r: Markup(svg_gauge(r))
    env.globals['status_badge'] = lambda s: Markup(status_badge(s))
    env.globals['highlight_seq'] = lambda seq, motif: Markup(highlight_motifs_in_sequence(seq, motif))

    template = env.get_template('longtr_pathogenic.html.jinja')

    return template.render(
        sample_name=sample_name,
        results=results,
        n_genotyped=sum(1 for r in results if r['genotyped']),
        n_pathogenic=sum(1 for r in results if r['locus_status'] == 'pathogenic'),
        n_intermediate=sum(1 for r in results if r['locus_status'] == 'intermediate'),
        n_uncertain=sum(1 for r in results if r['locus_status'] == 'uncertain'),
        n_normal=sum(1 for r in results if r['locus_status'] == 'normal'),
        n_missing=sum(1 for r in results if not r['genotyped']),
    )


def generate_report(vcf_path: str, strchive_json: str, longtr_bed: str, output: str):
    strchive = load_strchive_json(strchive_json)
    bed_entries = load_longtr_bed(longtr_bed)
    locus_index = build_locus_index(bed_entries)
    results, sample_name = scan_vcf(vcf_path, locus_index, strchive)
    html_content = generate_html(results, sample_name)
    with open(output, 'w') as f:
        f.write(html_content)
    def _fmt_ru(v):
        return f'{v:.0f}' if v == int(v) else f'{v:.1f}'

    print(f'Screened {len(results)} disease loci ({sum(1 for r in results if r["genotyped"])} genotyped)')
    for r in results:
        if r['locus_status'] in ('pathogenic', 'intermediate', 'uncertain'):
            print(f'  ⚠ {r["gene"]} ({r["disease"]}): {r["locus_status"]} — {_fmt_ru(r["allele1_ru"])}/{_fmt_ru(r["allele2_ru"])} repeats')
    print(f'Report written to {output}')


def cli_main():
    parser = ArgumentParser(
        description='Screen a LongTR VCF against STRchive disease-associated TR loci.',
    )
    parser.add_argument('--vcf_path', required=True, help='Path to LongTR VCF file')
    parser.add_argument(
        '--strchive_json',
        default=None,
        help='Path to STRchive-loci.json (downloaded automatically if not provided)',
    )
    parser.add_argument(
        '--longtr_bed',
        default=None,
        help='Path to STRchive LongTR BED catalog (downloaded automatically if not provided)',
    )
    parser.add_argument('--output', default='longtr_pathogenic.html', help='Output HTML file')
    args = parser.parse_args()

    cache_dir = Path.home() / '.cache' / 'longtr_pathogenic'
    cache_dir.mkdir(parents=True, exist_ok=True)

    json_path = args.strchive_json or str(download_if_missing(
        cache_dir / 'STRchive-loci.json', STRCHIVE_JSON_URL,
    ))
    bed_path = args.longtr_bed or str(download_if_missing(
        cache_dir / 'STRchive-disease-loci.hg38.longTR.bed', STRCHIVE_BED_URL,
    ))

    generate_report(args.vcf_path, json_path, bed_path, args.output)


if __name__ == '__main__':
    cli_main()
