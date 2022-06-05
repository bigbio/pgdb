/**
 * Add reference proteome
 */
process ADD_REFERENCE_PROTEOME {

    container "nfcore/pgdb:1.0.0"
    
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