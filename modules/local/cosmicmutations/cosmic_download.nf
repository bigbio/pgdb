/**
 * Download COSMIC Mutations
 */
process COSMIC_DOWNLOAD {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"
    
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
