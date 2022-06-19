/**
 * Download gnomAD variants (VCF) - requires gsutil
 */
process GNOMAD_DOWNLOAD {
    
    conda (params.enable_conda ? "conda-forge::gsutil=5.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"
        
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
