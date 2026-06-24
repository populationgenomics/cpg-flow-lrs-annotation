#!/usr/bin/env python3

"""
Check somalier self-relatedness results for a single participant.
If any SG pair has kinship below the threshold, sends a Slack alert.
Always exits 0.
"""

import csv
from argparse import ArgumentParser

from cpg_utils import config, slack
from loguru import logger


def run(
    pairs_fpath: str,
    sg_id: str,
    participant_id: str,
    dataset: str,
    kinship_threshold: float,
):
    logger.info(f'Checking self-relatedness for {sg_id} (participant {participant_id}) in {dataset}')
    logger.info(f'Kinship threshold: {kinship_threshold}')

    low_kinship_pairs = []

    with open(pairs_fpath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            relatedness = float(row['relatedness'])
            ibs0 = int(row['ibs0'])
            ibs2 = int(row['ibs2'])
            hom_concordance = float(row['hom_concordance'])
            hets_a = int(row['hets_a'])
            hets_b = int(row['hets_b'])
            shared_hets = int(row['shared_hets'])
            hom_alts_a = int(row['hom_alts_a'])
            hom_alts_b = int(row['hom_alts_b'])
            shared_hom_alts = int(row['shared_hom_alts'])
            n = int(row['n'])

            logger.info(
                f"  {row['#sample_a']} - {row['sample_b']}: "
                f'relatedness={relatedness}, ibs0={ibs0}, ibs2={ibs2}, '
                f'hom_concordance={hom_concordance}, '
                f'hets_a={hets_a}, hets_b={hets_b}, shared_hets={shared_hets}, '
                f'hom_alts_a={hom_alts_a}, hom_alts_b={hom_alts_b}, '
                f'shared_hom_alts={shared_hom_alts}, n={n}',
            )

            if relatedness < kinship_threshold:
                low_kinship_pairs.append({
                    'sample_a': row['#sample_a'],
                    'sample_b': row['sample_b'],
                    'relatedness': relatedness,
                    'ibs0': ibs0,
                    'ibs2': ibs2,
                })

    if not low_kinship_pairs:
        logger.info(
            f'{participant_id}: All pairs have kinship >= {kinship_threshold}',
        )
        return

    lines = [
        f'*[{dataset}] Self-relatedness check failed for {sg_id} (participant {participant_id})*',
        f'Expected kinship ~1.0 (threshold: {kinship_threshold}), found:',
    ]
    for pair in low_kinship_pairs:
        lines.append(
            f"  {pair['sample_a']} - {pair['sample_b']}: "
            f"kinship={pair['relatedness']}, "
            f"ibs0={pair['ibs0']}, ibs2={pair['ibs2']}",
        )

    text = '\n'.join(lines)
    logger.warning(text)

    if config.config_retrieve(
        ['workflow', 'somalier_self_check', 'send_to_slack'],
        default=True,
    ):
        slack.send_message(text)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--pairs-tsv', required=True)
    parser.add_argument('--sg-id', required=True)
    parser.add_argument('--participant-id', required=True)
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--kinship-threshold', type=float, default=0.9)
    args = parser.parse_args()
    run(
        pairs_fpath=args.pairs_tsv,
        sg_id=args.sg_id,
        participant_id=args.participant_id,
        dataset=args.dataset,
        kinship_threshold=args.kinship_threshold,
    )