/**
 * Creates the altORFs protein database
 */
process ADD_ALTORFS {

    container "nfcore/pgdb:1.0.0"
    
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