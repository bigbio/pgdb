/**
 * Generate proteindb from cosmic cell lines mutations
*/
process COSMIC_CELLLINES_PROTEINDB {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"

    when:
    params.cosmic_celllines

    input:
    file g
    file m
    file cosmic_config
    val cosmic_cellline_name

    output:
    path 'cosmic_celllines_proteinDB*.fa' ,emit:cosmic_celllines_proteindbs

    script:
    """
    pypgatk_cli.py cosmic-to-proteindb \\
        --config_file "$cosmic_config" \\
        --input_mutation $m \\
        --input_genes $g \\
        --filter_column 'Sample name' \\
        --accepted_values $cosmic_cellline_name \\
        --output_db cosmic_celllines_proteinDB.fa
    """
}
