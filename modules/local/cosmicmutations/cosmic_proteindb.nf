/**
 * Generate proteindb from cosmic mutations
*/
process COSMIC_PROTEINDB {
    
    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"

    when:
    params.cosmic

    input:
    file g
    file m
    file cosmic_config
    val cosmic_cancer_type

    output:
    path 'cosmic_proteinDB*.fa' ,emit: cosmic_proteindbs

    script:
    """
    pypgatk_cli.py cosmic-to-proteindb \\
        --config_file "$cosmic_config" \\
        --input_mutation $m --input_genes $g \\
        --filter_column 'Histology subtype 1' \\
        --accepted_values $cosmic_cancer_type \\
        --output_db cosmic_proteinDB.fa
    """
}
