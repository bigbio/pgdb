/**
 * Download gnomAD variants (VCF) - requires gsutil
 */
process GNOMAD_DOWNLOAD {
    
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
