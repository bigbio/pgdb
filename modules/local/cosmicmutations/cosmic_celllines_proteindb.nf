/**
 * Generate proteindb from cosmic cell lines mutations
*/
process COSMIC_CELLLINES_PROTEINDB {

    when:
    params.cosmic_celllines

    input:
    file g
    file m
    file cosmic_config

    output:
    path 'cosmic_celllines_proteinDB*.fa' ,emit:cosmic_celllines_proteindbs

    script:
    """
    pypgatk_cli.py cosmic-to-proteindb \\
        --config_file "$cosmic_config" \\
        --input_mutation $m \\
        --input_genes $g \\
        --filter_column 'Sample name' \\
        --accepted_values $params.cosmic_cellline_name \\
        --output_db cosmic_celllines_proteinDB.fa
    """
}