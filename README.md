# ![nf-core/pgdb](docs/images/nf-core-pgdb_logo.png)

The ProteoGenomics database generation workflow (**pgdb**) use the [pypgatk](https://github.com/bigbio/py-pgatk) and [nextflow](https://www.nextflow.io/) to create different protein databases for ProteoGenomics data analysis.

[![GitHub Actions CI Status](https://github.com/nf-core/pgdb/workflows/nf-core%20CI/badge.svg)](https://github.com/bigbio/pgdb/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/pgdb/workflows/nf-core%20linting/badge.svg)](https://github.com/bigbio/pgdb/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/pgdb.svg)](https://hub.docker.com/r/nfcore/pgdb)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23pgdb-4A154B?logo=slack)](https://nfcore.slack.com/channels/pgdb)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**nf-core/pgdb** is a bioinformatics best-practise analysis pipeline for

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/pgdb -profile test,<docker/singularity/podman/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run nf-core/pgdb -profile <docker/singularity/podman/conda/institute> --ensembl_name homo_sapines --ensembl false
    ```

See [usage docs](https://nf-co.re/pgdb/usage) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following:

![ProteoGenomics Database](/docs/images/pgdb-databases.png)

* Download protein databases from ENSEMBL
* Translate from Genomics Variant databases into ProteoGenomics Databases (`COSMIC`, `GNOMAD`)
* Add to a Reference proteomics database, non-coding RNAs + pseudogenes.
* Compute Decoy for a proteogenomics databases

## Documentation

The nf-core/pgdb pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/pgdb/usage) and [output](https://nf-co.re/pgdb/output).

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

## Credits

nf-core/pgdb was originally written by Husen M. Umer & Yasset Perez-Riverol.

We thank the following people for their extensive assistance in the development
of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#pgdb` channel](https://nfcore.slack.com/channels/pgdb) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/pgdb for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

In addition, references of tools and data used in this pipeline are as follows:

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
