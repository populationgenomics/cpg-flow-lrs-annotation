# CPG-Flow - Long Read Sequencing Annotation for seqr

Version 0.2.1

A CPG workflow for creating annotated callsets from long read data, using the [cpg-flow](https://github.com/populationgenomics/cpg-flow) pipeline framework.

## Purpose

This workflow is designed to process long-read sequencing data and create callsets compatible with [seqr](https://github.com/populationgenomics/seqr). It automates the steps required to query, reformat, annotate, and export VCF files derived from long-read sequencing. It also supports the conversion of BAM files to CRAMs.

It is intended to be used for SNPs/Indels VCFs and SVs VCFs from different sequencing technologies, such as PacBio, Oxford Nanopore, and others. The inputs can be configured by the user, allowing for flexibility in the types of long-read data queried and processed.

## Workflow Overview

### Annotation

1. Query Metamist for long-read sequencing groups and their VCF analyses based on the filters specified in the configuration
   2. If the SGs from the input cohorts do not have VCFs matching the filter criteria, the workflow will fail.
2. Perform necessary reformatting, reheadering, and normalization of the VCFs
3. Merge the VCFs and annotate the merged VCF with VEP (for SNPs Indels) or STRVCTVRE (for SVs)
4. Write the annotated VCF to a Matrix table
5. Export the Matrix table to an elasticsearch index

### Conversion

1. Query Metamist for long-read sequencing groups and their BAM assays
2. Convert BAM files to CRAM files using `samtools`

## Directory Structure

```commandline
src
├── lrs_annotation
│   ├── __init__.py
│   ├── run_workflow.py
│   ├── lrs_annotation.toml
│   ├── bam_to_cram_stages.py
│   ├── bam_to_cram_stages.toml
│   ├── snps_indels_annotation_stages.py
│   ├── snps_indels_annotation.toml
│   ├── svs_annotation_stages.py
│   ├── svs_annotation.toml
│   ├── inputs.py
│   ├── utils.py
│   ├── jobs
│   │   ├── snps_indels
│   │   │   ├── AnnotateCohortMatrixtable.py
│   │   │   └── AnnotateDatasetMatrixtable.py
│   │   │   └── ...
│   │   └── svs
│   │   │   ├── AnnotateCohortMatrixtable.py
│   │   │   └── AnnotateDatasetMatrixtable.py
│   │   │   └── ...
│   ├── scripts
│   │   ├── snps_indels
│   │   │   ├── __init__.py
│   │   │   ├── annotate_cohort_snps_indels.py
│   │   │   ├── ...
│   │   ├── svs
│   │   │   ├── __init__.py
│   │   │   ├── annotate_cohort_svs.py
│   │   │   ├── ...
│   ├── hail_scripts/computed_fields
│   │   ├── __init__.py
│   │   ├── variant_id.py
│   │   └── vep.py
```

## Key Components

`lrs_annotation.toml` contains the main configuration for the workflow, including a number of mandatory options shared between the SNPs/Indels and SVs workflows.

`snps_indels_annotation_stages.py` and `svs_annotation_stages.py` contain Stages for the workflows, with the actual logic imported from files in `jobs`.

`snps_indels_annotation.toml` and `svs_annotation.toml` are config files to be submitted with the workflow via analysis-runner. They contain settings which are used to configure the workflow, such as filters for fetching input VCFs, resource allocation for specific jobs, and references.

`scripts/` contains the scripts that are required for the workflows, often as Query on Batch jobs.

`hail_scripts/computed_fields/` contains computed fields required for formatting the callsets for seqr, borrowed from the [seqr-loading-pipelines](https://github.com/broadinstitute/seqr-loading-pipelines/tree/master/hail_scripts/computed_fields)

`inputs.py` contains functions to query Metamist for long-read sequencing groups and their VCF analyses, as well as functions to fetch the necessary VCFs based on the configuration.

`utils.py` contains utility functions used across the workflows, such as parsing command line arguments, submitting cromwell jobs, and reading the configuration.
