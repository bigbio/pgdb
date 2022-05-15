/**
 * Generate proteindb from local cosmic cell lines mutations
*/
process COSMIC_CELLLINES_PROTEINDB_LOCAL {

    when:
    params.cosmiccelllines_genes && params.cosmiccelllines_mutations

    input:
    if (params.cosmiccelllines_genes&&params.cosmiccelllines_mutations) {
        file cosmiccelllines_genes
        file cosmiccelllines_mutations
    }
    
    file cosmic_config

    output:
    path 'cosmic_celllines_proteinDB*.fa' ,emit:cosmic_celllines_proteindbs_uselocal

    script:
    """
    pypgatk_cli.py cosmic-to-proteindb \\
        --config_file "$cosmic_config" \\
        --input_mutation $cosmiccelllines_mutations \\
        --input_genes $cosmiccelllines_genes \\
        --filter_column 'Sample name' \\
        --accepted_values $params.cosmic_cellline_name \\
        --output_db cosmic_celllines_proteinDB.fa
    """
}