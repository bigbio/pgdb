/**
 * Download COSMIC Cell Lines Mutations
 */
process COSMIC_CELLLINES_DOWNLOAD {

    when:
    params.cosmic_celllines

    output:
    path "database_cosmic/All_CellLines_Genes.fasta" ,emit: cosmic_celllines_genes
    path "database_cosmic/CosmicCLP_MutantExport.tsv" ,emit:cosmic_celllines_mutations

    script:
    base64 = ''
    """
    base64=`echo "$params.cosmic_user_name:$params.cosmic_password" | base64`
    curl -o database_cosmic/All_CellLines_Genes.fasta.gz --create-dirs `curl -H "Authorization: \$base64" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cell_lines/$params.cosmic_dababase_version/All_CellLines_Genes.fasta.gz |grep -Po 'url[" :]+\\K[^"]+'`
    curl -o database_cosmic/CosmicCLP_MutantExport.tsv.gz --create-dirs `curl -H "Authorization: \$base64" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cell_lines/$params.cosmic_dababase_version/CosmicCLP_MutantExport.tsv.gz |grep -Po 'url[" :]+\\K[^"]+'`
    gunzip database_cosmic/*.gz
    """
}
