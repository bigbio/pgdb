/**
 * Download VCF files from ensembl for the particular species.
 */
process ENSEMBL_VCF_DOWNLOAD {

    when:
    params.ensembl

    input:
    file ensembl_downloader_config

    output:
    path "database_ensembl/*.vcf" ,emit: ensembl_vcf_files

    script:
    """
    pypgatk_cli.py ensembl-downloader \\
        --config_file $ensembl_downloader_config \\
        --ensembl_name $params.ensembl_name \\
        -sg -sp -sc -sd -sn
    """
}
