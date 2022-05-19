/**
 * Generate protein databse for a given VCF
 */
process GTF_TO_FASTA {

    when:
    params.vcf

    input:
    file g
    file f

    output:
    path "transcripts.fa" ,emit:gtf_transcripts_fasta

    script:
    """
    gffread -w transcripts.fa -g $f $g
    """
}