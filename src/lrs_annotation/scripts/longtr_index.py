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
                family=item.get('family', ''),
                ext_participant=item.get('ext_participant', ''),
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


def main(manifest: str, dataset_name: str, output: str) -> None:
    report_items, loci_lists = load_manifest(manifest)
    entries = build_entries_from_reports(report_items)

    template_dir = Path(__file__).resolve().parent / 'templates'
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(template_dir)),
        autoescape=True,
    )
    template = env.get_template('longtr_index.html.jinja')
    content = template.render(reports=entries, dataset=dataset_name, loci_lists=loci_lists)

    with Path(output).open('w') as f:
        f.write(content)


if __name__ == '__main__':
    parser = ArgumentParser(description='Generate an index page for LongTR pathogenic reports')
    parser.add_argument('--manifest', required=True, help='JSON manifest listing all reports')
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument('--output', required=True, help='Output HTML file path')
    args = parser.parse_args()
    main(manifest=args.manifest, dataset_name=args.dataset, output=args.output)
