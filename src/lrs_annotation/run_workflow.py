#!/usr/bin/env python3

"""
Choose a workflow to run for LRS annotation.
"""

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow

from lrs_annotation.bam_to_cram_stages import ConvertBamToCram
from lrs_annotation.snps_indels_annotation_stages import ExportSnpsIndelsMtToESIndex
from lrs_annotation.svs_annotation_stages import ExportSVsMtToElasticIndex


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    parser.add_argument('--workflow', type=str, choices=['bam_to_cram', 'snps_indels', 'svs'], default='snps_indels')
    args = parser.parse_args()

    if args.workflow == 'bam_to_cram':
        stages = [ConvertBamToCram]
    elif args.workflow == 'snps_indels':
        stages = [ExportSnpsIndelsMtToESIndex]
    elif args.workflow == 'svs':
        stages = [ExportSVsMtToElasticIndex]
    else:
        raise ValueError(f'Unknown workflow: {args.workflow}')

    run_workflow(name='lrs_annotation', stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
