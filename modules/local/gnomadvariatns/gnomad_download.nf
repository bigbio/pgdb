/**
 * Download gnomAD variants (VCF) - requires gsutil
 */
process GNOMAD_DOWNLOAD {

    container "nfcore/pgdb:1.0.0"
    
    when:
    params.gnomad

    input:
    val g

    output:
    path "*.vcf.bgz" ,emit:gnomad_vcf_bgz

    script:
    """
    gsutil cp $g .
    """
}
