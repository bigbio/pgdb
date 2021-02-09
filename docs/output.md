# nf-core/pgdb: Output

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/pgdb/output](https://nf-co.re/pgdb/output)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This document describes the output produced by the pipeline. The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Ensembl](#ensembl) - Download the Ensembl databases proteins, GTF, etc
* [Variants](#variants) - Add the COSMIC and cPortal variant databases
* [Output](#output) - Output results including clean databases and decoy generation

## Ensembl

The pipeline will download the the ENSEMBL protein reference proteome, this will be added to the result files.

## Variants

The COSMIC or cBioPortal variants.

## Output

 /fasta_database.fa

 The FASTA database including all the protein sequences including the reference proteomes, variants, pseudo-genes, etc.
