/**
 * Download COSMIC Mutations
 */
process COSMIC_DOWNLOAD {

    when:
    params.cosmic || params.cosmic_celllines

    input:
    file cosmic_config

    output:
    path "database_cosmic/All_COSMIC_Genes.fasta" ,emit:cosmic_genes
    path "database_cosmic/CosmicMutantExport.tsv" ,emit:cosmic_mutations
    path "database_cosmic/All_CellLines_Genes.fasta" ,emit:cosmic_celllines_genes
    path "database_cosmic/CosmicCLP_MutantExport.tsv" ,emit:cosmic_celllines_mutations

    script:
    """
    pypgatk_cli.py cosmic-downloader \\
        --config_file "$cosmic_config" \\
        --username $params.cosmic_user_name \\
        --password $params.cosmic_password
    """
}