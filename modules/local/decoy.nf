/**
 * Create the decoy database using DecoyPYrat
 * Decoy sequences will have "DECOY_" prefix tag to the protein accession.
 */
process DECOY {

    container "nfcore/pgdb:1.0.0"
    
    publishDir "${params.outdir}/", mode: params.publish_dir_mode,
        saveAs: { filename -> params.final_database_protein }

    when:
    params.decoy

    input:
    file f
    file protein_decoy_config

    output:
    path 'decoy_database.fa' ,emit: fasta_decoy_db_ch

    script:
    """
    pypgatk_cli.py generate-decoy \\
        --method "$params.decoy_method" \\
        --enzyme "$params.decoy_enzyme" \\
        --config_file $protein_decoy_config \\
        --input_database $f \\
        --decoy_prefix "$params.decoy_prefix" \\
        --output_database decoy_database.fa
    """
}
