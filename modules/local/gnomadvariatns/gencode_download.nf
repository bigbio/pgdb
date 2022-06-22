/**
 * Download gencode files (fasta and gtf)
 */
process GENCODE_DOWNLOAD {
    
    when:
    params.gnomad

    input:
    val g

    output:
    path("gencode.v19.pc_transcripts.fa") ,emit:gencode_fasta
    path("gencode.v19.annotation.gtf") ,emit:gencode_gtf

    script:
    """
    wget ${g}/gencode.v19.pc_transcripts.fa.gz
    wget ${g}/gencode.v19.annotation.gtf.gz
    gunzip *.gz
    """
}
