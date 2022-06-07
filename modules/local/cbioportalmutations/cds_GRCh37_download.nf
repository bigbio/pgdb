/**
 * Download GRCh37 CDS file from ENSEMBL release 75
 */
process CDS_GRCH37_DOWNLOAD {

    
    when:
    params.cbioportal

    output:
    path("Homo_sapiens.GRCh37.75.cds.all.fa") ,emit:ch_GRCh37_cds

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.75.cds.all.fa.gz
    gunzip *.gz
    """
}

