/**
 * Concatenate cDNA and ncRNA databases
 **/
process MERGE_CDNAS {
    
    input:
    file a
    file b

    output:
    path 'total_cdnas.fa',emit: total_cdnas

    script:
    """
    cat $a >> total_cdnas.fa
    cat $b >> total_cdnas.fa
    """
}
