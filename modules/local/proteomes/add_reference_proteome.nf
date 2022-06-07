/**
 * Add reference proteome
 */
process ADD_REFERENCE_PROTEOME {
    
    when:
    params.add_reference

    input:
    file reference_proteome

    output:
    path 'reference_proteome.fa', emit: ensembl_protein_database

    script:
    """
    cat $reference_proteome >> reference_proteome.fa
    """
}
