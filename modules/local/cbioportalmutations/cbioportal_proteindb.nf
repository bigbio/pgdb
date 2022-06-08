/**
 * Generate proteinDB from cBioPortal mutations
 */
process CBIOPORTAL_PROTEINDB {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"

    when:
    params.cbioportal

    input:
    file g
    file m
    file s
    file cbioportal_config
    val cbioportal_filter_column
    val cbioportal_accepted_values

    output:
    path 'cbioPortal_proteinDB*.fa' ,emit: cBioportal_proteindb

    script:
    """
    pypgatk_cli.py cbioportal-to-proteindb \\
        --config_file $cbioportal_config \\
        --input_mutation $m \\
        --input_cds $g \\
        --clinical_sample_file $s \\
        --filter_column $cbioportal_filter_column \\
        --accepted_values $cbioportal_accepted_values \\
        --output_db cbioPortal_proteinDB.fa
    """
}
