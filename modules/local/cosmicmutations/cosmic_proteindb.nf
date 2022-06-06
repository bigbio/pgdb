/**
 * Generate proteindb from cosmic mutations
*/
process COSMIC_PROTEINDB {
    
    container "nfcore/pgdb:1.0.0"

    when:
    params.cosmic

    input:
    file g
    file m
    file cosmic_config

    output:
    path 'cosmic_proteinDB*.fa' ,emit: cosmic_proteindbs

    script:
    """
    pypgatk_cli.py cosmic-to-proteindb \\
        --config_file "$cosmic_config" \\
        --input_mutation $m --input_genes $g \\
        --filter_column 'Histology subtype 1' \\
        --accepted_values $params.cosmic_cancer_type \\
        --output_db cosmic_proteinDB.fa
    """
}
