/**
 * Generate proteinDB from cBioPortal mutations
 */
process CBIOPORTAL_PROTEINDB {

    container "nfcore/pgdb:1.0.0"

    when:
    params.cbioportal

    input:
    file g
    file m
    file s
    file cbioportal_config

    output:
    path 'cbioPortal_proteinDB*.fa' ,emit: cBioportal_proteindb

    script:
    """
    pypgatk_cli.py cbioportal-to-proteindb \\
        --config_file $cbioportal_config \\
        --input_mutation $m \\
        --input_cds $g \\
        --clinical_sample_file $s \\
        --filter_column $params.cbioportal_filter_column \\
        --accepted_values $params.cbioportal_accepted_values \\
        --output_db cbioPortal_proteinDB.fa
    """
}
