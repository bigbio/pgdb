/**
 * Concatenate all generated databases from merged_databases channel to the final_database_protein file
 */
process MERGE_PROTEINDBS {
    
    publishDir "${params.outdir}/", mode: params.publish_dir_mode,
        // Final step if not cleaning or creating a decoy database - save output to params.final_database_protein
        saveAs: { filename ->
            params.clean_database || params.decoy ? null : params.final_database_protein
        }

    input:
    file("proteindb*")

    output:
    path 'merged_databases.fa' ,emit: to_clean_ch

    script:
    """
    cat proteindb* > merged_databases.fa
    """
}
