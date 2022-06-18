/*
 * Parse software version numbers
 */
process GETSOFTWAREVERSIONS {

    output:
    path "v_pipeline.txt" ,emit: v_pipeline
    path "v_nextflow.txt" ,emit: v_nextflow
    path "pypgatk_version.txt" ,emit: pypgatk_version

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    pip list |grep pypgatk > pypgatk_version.txt
    """
}
