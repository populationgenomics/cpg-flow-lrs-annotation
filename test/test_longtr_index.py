"""Tests for the dataset-level LongTR index page builder."""

import json

import pytest

from lrs_annotation.scripts.longtr_index import (
    build_entries_from_reports,
    enrich_manifest_from_json,
    parse_json_entries,
)


def test_parse_json_entries_splits_on_first_two_colons_only():
    """Triples key on (sg_id, report_type); any further colons stay part of the path."""
    entries = [
        'CPG1:neuro:/io/batch/abc/a.json',
        'CPG2:paediatric:/io/batch/weird:name/b.json',
    ]
    assert parse_json_entries(entries) == {
        ('CPG1', 'neuro'): '/io/batch/abc/a.json',
        ('CPG2', 'paediatric'): '/io/batch/weird:name/b.json',
    }


@pytest.fixture
def report_json(tmp_path):
    """A per-SG JSON report containing one locus of each flagged status plus a normal one."""
    path = tmp_path / 'report.json'
    path.write_text(
        json.dumps(
            {
                'sample_name': 'S1',
                'summary': {},
                'loci': [
                    {'gene': 'TCF4', 'locus_status': 'pathogenic', 'genotyped': True},
                    {'gene': 'RFC1', 'locus_status': 'intermediate', 'genotyped': True},
                    {'gene': 'POLG', 'locus_status': 'uncertain', 'genotyped': True},
                    {'gene': 'HTT', 'locus_status': 'normal', 'genotyped': True},
                    {'gene': 'FMR1', 'locus_status': 'not_genotyped', 'genotyped': False},
                ],
            }
        )
    )
    return str(path)


def test_enrich_manifest_extracts_flagged_and_missing(report_json):
    """Every flagged status reaches the index; the filter is derived from STATUS_COLORS."""
    items = [{'sample': 'CPG1', 'report_type': 'neuro'}]
    enrich_manifest_from_json(items, {('CPG1', 'neuro'): report_json})

    flagged = {f['gene']: f['status'] for f in items[0]['flagged_loci']}
    assert flagged == {
        'TCF4': 'pathogenic',
        'RFC1': 'intermediate',
        'POLG': 'uncertain',
    }
    # normal loci are not flagged, and ungenotyped loci are reported separately
    assert 'HTT' not in flagged
    assert items[0]['missing_loci'] == ['FMR1']


def test_build_entries_groups_loci_by_colour(report_json):
    """Flagged loci are grouped by the colour the index template renders them with."""
    items = [{'sample': 'CPG1', 'report_type': 'neuro'}]
    enrich_manifest_from_json(items, {('CPG1', 'neuro'): report_json})
    entry = build_entries_from_reports(items)[0]

    assert entry.loci_of_interest == {
        'Red': ['TCF4'],
        'Orange': ['RFC1'],
        'Grey': ['POLG'],
    }


def test_enrich_manifest_tolerates_missing_json(tmp_path):
    """A missing JSON report leaves the entry unenriched rather than aborting the index build."""
    items = [{'sample': 'CPG1', 'report_type': 'neuro'}]
    enrich_manifest_from_json(items, {('CPG1', 'neuro'): str(tmp_path / 'absent.json')})
    assert 'flagged_loci' not in items[0]
