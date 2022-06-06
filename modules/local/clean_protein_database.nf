/**
 * clean the database for stop codons, and unwanted AA like: *, also remove proteins with less than 6 AA
 */
process CLEAN_PROTEIN_DATABASE {

    container "nfcore/pgdb:1.0.0"

    publishDir "${params.outdir}/", mode: params.publish_dir_mode,
        // Final step if not creating a decoy database - save output to params.final_database_protein
        saveAs: { filename ->
            params.decoy ? null : params.final_database_protein
        }

    when:
    params.clean_database

    input:
    file file
    file e

    output:
    path 'database_clean.fa' ,emit: clean_database_sh

    script:
    stop_codons = ''
    if (params.add_stop_codons){
       stop_codons = "--add_stop_codons"
    }

    """
    pypgatk_cli.py ensembl-check \\
        -in "$file" \\
        --config_file "$e" \\
        -out database_clean.fa \\
        --num_aa "$params.minimum_aa" \\
        "$stop_codons"
    """
}
