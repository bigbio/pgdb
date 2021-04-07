# nf-core/pgdb: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - Lovely Logan [18.10.2020]

Initial release of nf-core/pgdb, created with the [nf-core](https://nf-co.re/) template.

### `Added`

The initial version of the pipeline features the following steps:

    - (optional) ENSEMBL Reference proteomes included in final proteome
    - (optional) Generate the decoy database and attach it to the final proteome
    - Convert a Variant genome database like COSMIC or CBioPortal to proteomes
    -

### `Known issues`

If you experience nextflow running forever after a failed step, try settings errorStrategy = terminate. See the corresponding [nextflow issue](https://github.com/nextflow-io/nextflow/issues/1457).

### `Fixed`

### `Dependencies`

### `Deprecated`
