/**
 * Generate proteindb from local cosmic mutations
*/
process COSMIC_PROTEINDB_LOCAL {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"
    
    when:
    params.cosmicgenes && params.cosmicmutations

    input:
    file cosmic_config
    val cosmic_cancer_type
    file cosmicmutations
    file cosmicgenes

    output:
    file 'cosmic_proteinDB*.fa' into cosmic_proteindbs_uselocal

    script:
    """
    pypgatk_cli.py cosmic-to-proteindb \\
        --config_file "$cosmic_config" \\
        --input_mutation $cosmicmutations --input_genes $cosmicgenes \\
        --filter_column 'Histology subtype 1' \\
        --accepted_values $cosmic_cancer_type \\
        --output_db cosmic_proteinDB.fa
    """
}
