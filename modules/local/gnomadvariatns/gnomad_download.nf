/**
 * Download gnomAD variants (VCF) - requires gsutil
 */
process GNOMAD_DOWNLOAD {
    
    conda (params.enable_conda ? "conda-forge::gsutil=5.10" : null)

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
