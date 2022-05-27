/**
 * Generate proteindb from local cosmic mutations
*/
process COSMIC_PROTEINDB_LOCAL {

    when:
    params.cosmicgenes && params.cosmicmutations

    input:
    if (params.cosmicgenes&&params.cosmicmutations) {
        file cosmicgenes
        file cosmicmutations
    }
    file cosmic_config

    output:
    file 'cosmic_proteinDB*.fa' into cosmic_proteindbs_uselocal

    script:
    """
    pypgatk_cli.py cosmic-to-proteindb \\
        --config_file "$cosmic_config" \\
        --input_mutation $cosmicmutations --input_genes $cosmicgenes \\
        --filter_column 'Histology subtype 1' \\
        --accepted_values $params.cosmic_cancer_type \\
        --output_db cosmic_proteinDB.fa
    """
}