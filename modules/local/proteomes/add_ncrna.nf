/**
 * Creates the ncRNA protein database
 */
process ADD_NCRNA {

    container "nfcore/pgdb:1.0.0"
    
    when:
    params.ncrna

    input:
    file x
    file ensembl_config

    output:
    path 'ncRNAs_proteinDB.fa',emit: optional_ncrna

    script:
    """
    pypgatk_cli.py dnaseq-to-proteindb \\
        --config_file "$ensembl_config" \\
        --input_fasta $x \\
        --output_proteindb ncRNAs_proteinDB.fa \\
        --include_biotypes "${params.biotypes['ncRNA']}" \\
        --skip_including_all_cds --var_prefix ncRNA_
    """
}