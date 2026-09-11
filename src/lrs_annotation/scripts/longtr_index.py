"""
Generate a dataset-level index page linking to individual LongTR pathogenic
TR expansion reports, with filtering and loci-of-interest display.
"""

import json
import re
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import jinja2


@dataclass
class IndexEntry:
    sample: str
    family: str
    ext_participant: str
    ext_sample: str
    affected_status: str
    report_type: str
    run_date: str
    missing_loci: str
    url: str
    loci_of_interest: dict[str, list[str]] = field(default_factory=dict)

    def __key(self) -> tuple[str, str, str]:
        return self.sample, self.report_type, self.url

    def __hash__(self) -> int:
        return hash(self.__key())


STATUS_COLORS = {
    'pathogenic': 'Red',
    'intermediate': 'Orange',
    'uncertain': 'Grey',
}


def load_manifest(manifest_path: str) -> tuple[list[dict], dict[str, list[str]]]:
    with open(manifest_path) as f:
        raw = json.load(f)

    if isinstance(raw, list):
        return raw, {}
    return raw.get('reports', []), raw.get('loci_lists', {})


def load_json_map(json_map_path: str) -> dict[tuple[str, str], str]:
    """Parse a TSV mapping of sg_id, report_type, json_path."""
    mapping: dict[tuple[str, str], str] = {}
    with open(json_map_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            min_tsv_columns = 3
            if len(parts) >= min_tsv_columns:
                mapping[(parts[0], parts[1])] = parts[2]
    return mapping


def enrich_manifest_from_json(
    report_items: list[dict],
    json_map: dict[tuple[str, str], str],
) -> None:
    """Read JSON report files and add flagged_loci/missing_loci to manifest entries in-place."""
    for item in report_items:
        key = (item['sample'], item.get('report_type', 'default'))
        json_path = json_map.get(key)
        if not json_path:
            continue
        try:
            with open(json_path) as f:
                report = json.load(f)
            flagged = []
            missing = []
            for locus in report.get('loci', []):
                status = locus.get('locus_status', 'normal')
                if status in ('pathogenic', 'intermediate', 'uncertain'):
                    flagged.append({'gene': locus['gene'], 'status': status})
                if not locus.get('genotyped', True):
                    missing.append(locus['gene'])
            item['flagged_loci'] = flagged
            item['missing_loci'] = missing
        except (OSError, json.JSONDecodeError, KeyError) as e:
            print(f'Warning: could not read {json_path}: {e}')


def build_entries_from_reports(report_items: list[dict]) -> list[IndexEntry]:
    entries = []
    for item in report_items:
        loci_of_interest: dict[str, list[str]] = defaultdict(list)
        for locus in item.get('flagged_loci', []):
            color = STATUS_COLORS.get(locus.get('status', ''), 'Grey')
            loci_of_interest[color].append(locus['gene'])

        entries.append(
            IndexEntry(
                sample=item['sample'],
                family=item.get('family_id', item.get('family', '')),
                ext_participant=item.get('external_id', item.get('ext_participant', '')),
                ext_sample=item.get('ext_sample', ''),
                affected_status=item.get('affected_status', ''),
                report_type=re.sub(r'[-_]', ' ', item.get('report_type', 'default')).title(),
                run_date=item.get('date', ''),
                missing_loci=', '.join(item.get('missing_loci', [])),
                url=item.get('url', ''),
                loci_of_interest=dict(loci_of_interest),
            ),
        )
    return entries


def main(manifest: str, dataset_name: str, output: str, json_map_path: str | None = None) -> None:
    report_items, loci_lists = load_manifest(manifest)
    if json_map_path:
        enrich_manifest_from_json(report_items, load_json_map(json_map_path))
    entries = build_entries_from_reports(report_items)

    template_dir = Path(__file__).resolve().parent / 'templates'
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(template_dir)),
        autoescape=True,
    )
    env.filters['gene_name'] = lambda locus_id: locus_id.rsplit('_', 1)[-1] if '_' in locus_id else locus_id
    template = env.get_template('longtr_index.html.jinja')
    content = template.render(reports=entries, dataset=dataset_name, loci_lists=loci_lists)

    with Path(output).open('w') as f:
        f.write(content)


if __name__ == '__main__':
    parser = ArgumentParser(description='Generate an index page for LongTR pathogenic reports')
    parser.add_argument('--manifest', required=True, help='JSON manifest listing all reports')
    parser.add_argument(
        '--json-map', dest='json_map', default=None, help='TSV mapping sg_id/report_type to JSON report paths'
    )
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument('--output', required=True, help='Output HTML file path')
    args = parser.parse_args()
    main(manifest=args.manifest, dataset_name=args.dataset, output=args.output, json_map_path=args.json_map)
