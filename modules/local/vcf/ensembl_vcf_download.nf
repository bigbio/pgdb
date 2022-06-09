/**
 * Download VCF files from ensembl for the particular species.
 */
process ENSEMBL_VCF_DOWNLOAD {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"
    
    when:
    params.ensembl

    input:
    file ensembl_downloader_config
    val ensembl_name

    output:
    path "database_ensembl/*.vcf" ,emit: ensembl_vcf_files

    script:
    """
    pypgatk_cli.py ensembl-downloader \\
        --config_file $ensembl_downloader_config \\
        --ensembl_name $ensembl_name \\
        -sg -sp -sc -sd -sn
    """
}

