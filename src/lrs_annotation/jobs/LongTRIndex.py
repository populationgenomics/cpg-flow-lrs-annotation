"""
Job to generate an index HTML page linking to all LongTR pathogenic reports
for a dataset, with sample metadata from metamist.
"""

import json
import re

import loguru

from hailtop.batch.job import Job

from cpg_utils import Path, config, hail_batch
from metamist.graphql import gql, query

METADATA_QUERY = gql(
    """
    query SgMetadata($project: String!, $sgIds: [String!]!) {
        project(name: $project) {
            meta
            sequencingGroups(id: {in_: $sgIds}) {
                id
                sample {
                    participant {
                        externalId
                        families {
                            externalId
                        }
                        familyParticipants {
                            affected
                        }
                    }
                }
            }
        }
    }
    """
)


def get_sg_metadata(
    dataset: str,
    sg_ids: list[str],
) -> tuple[dict[str, dict[str, str | int]], str]:
    """
    Query metamist for per-SG metadata (family ID, external ID, affected status)
    and the project display name.
    """
    query_dataset = dataset
    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    result = query(METADATA_QUERY, variables={'project': query_dataset, 'sgIds': sg_ids})

    project = result.get('project', {})
    display_name = project.get('meta', {}).get('display_name', dataset)

    sg_metadata: dict[str, dict[str, str | int]] = {}
    for group in project.get('sequencingGroups', []):
        sg_id = group.get('id')
        try:
            participant = group['sample']['participant']
            sg_metadata[sg_id] = {
                'family_id': participant['families'][0]['externalId'],
                'external_id': participant['externalId'],
                'affected': participant['familyParticipants'][0]['affected'],
            }
        except (KeyError, IndexError, TypeError):
            if sg_id in sg_ids:
                loguru.logger.warning(f'Missing metadata for {sg_id}')
            continue

    return sg_metadata, display_name


def _affected_label(affected: int | str) -> str:
    match affected:
        case 1:
            return 'Unaffected'
        case 2:
            return 'Affected'
        case _:
            return 'Unknown'


def longtr_index_page(
    dataset_name: str,
    sg_report_outputs: dict[str, dict[str, Path]],
    loci_lists: dict[str, list[str]],
    output_path: Path,
    job_attrs: dict[str, str],
) -> Job:
    """
    Generate an index HTML page linking to all LongTR pathogenic reports for a dataset.
    Queries metamist for sample metadata (family, external ID, affected status).
    """
    from lrs_annotation.scripts import longtr_index  # noqa: PLC0415

    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_job(f'LongTR Index Page for {dataset_name}', job_attrs | {'tool': 'longtr_index'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    sg_ids = list(sg_report_outputs.keys())
    sg_metadata, display_name = get_sg_metadata(dataset_name, sg_ids)
    dataset_title = re.sub(r'[-_]', ' ', display_name).title()

    file_prefix = config.config_retrieve(['storage', dataset_name, 'web'])
    html_prefix = config.config_retrieve(['storage', dataset_name, 'web_url'])

    manifest_data = []
    for sg_id, output_dict in sg_report_outputs.items():
        meta = sg_metadata.get(sg_id, {})
        family_id = meta.get('family_id', '')
        external_id = meta.get('external_id', '')
        affected_status = _affected_label(meta.get('affected', 0))

        for key, report_path in output_dict.items():
            if not key.endswith('_html') and key != 'html':
                continue
            report_type = key.removesuffix('_html') if key != 'html' else 'default'
            url = str(report_path).replace(file_prefix, html_prefix)
            manifest_data.append(
                {
                    'sample': sg_id,
                    'family_id': family_id,
                    'external_id': external_id,
                    'affected_status': affected_status,
                    'report_type': report_type,
                    'url': url,
                }
            )

    manifest = {'reports': manifest_data, 'loci_lists': loci_lists}
    job.command(f"""
    cat > {job.manifest} << 'MANIFEST_EOF'
{json.dumps(manifest, indent=2)}
MANIFEST_EOF

    python3 {longtr_index.__file__} \\
        --manifest {job.manifest} \\
        --dataset '{dataset_title}' \\
        --output {job.output}
    """)

    batch_instance.write_output(job.output, str(output_path))
    return job
