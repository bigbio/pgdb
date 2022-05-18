/**
 * Creates the pseudogenes protein database
 */
process ADD_PSEUDOGENES {

    when:
    params.pseudogenes

    input:
    file x
    file ensembl_config

    output:
    path 'pseudogenes_proteinDB.fa' ,emit: optional_pseudogenes

    script:
    """
    pypgatk_cli.py dnaseq-to-proteindb \\
        --config_file "$ensembl_config" \\
        --input_fasta "$x" \\
        --output_proteindb pseudogenes_proteinDB.fa \\
        --include_biotypes "${params.biotypes['pseudogene']}" \\
        --skip_including_all_cds \\
        --var_prefix pseudo_
    """
}
