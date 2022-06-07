/*
 * Parse software version numbers
 */
process GET_SOFTWARE_VERSIONS {
    
    publishDir "${params.outdir}/pipeline_info", 
        mode: params.publish_dir_mode,
        saveAs: { filename ->  if (filename.indexOf('.csv') > 0) filename else null }

    output:
    path 'software_versions_mqc.yaml' ,emit: ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}
