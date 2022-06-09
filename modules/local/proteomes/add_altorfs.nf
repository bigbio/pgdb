/**
 * Creates the altORFs protein database
 */
process ADD_ALTORFS {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"
    
    when:
    params.altorfs

    input:
    file x
    file ensembl_config

    output:
    path 'altorfs_proteinDB.fa' ,emit:optional_altorfs

    script:
    """
    pypgatk_cli.py dnaseq-to-proteindb \\
        --config_file "$ensembl_config" \\
        --input_fasta "$x" \\
        --output_proteindb altorfs_proteinDB.fa \\
        --include_biotypes "${params.biotypes['protein_coding']}'" \\
        --skip_including_all_cds \\
        --var_prefix altorf_
    """
}
