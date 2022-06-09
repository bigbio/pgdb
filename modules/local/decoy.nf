/**
 * Create the decoy database using DecoyPYrat
 * Decoy sequences will have "DECOY_" prefix tag to the protein accession.
 */
process DECOY {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"
    
    publishDir "${params.outdir}/", mode: params.publish_dir_mode,
        saveAs: { filename -> params.final_database_protein }

    when:
    params.decoy

    input:
    file f
    file protein_decoy_config
    val decoy_method
    val decoy_enzyme
    val decoy_prefix

    output:
    path 'decoy_database.fa' ,emit: fasta_decoy_db_ch

    script:
    """
    pypgatk_cli.py generate-decoy \\
        --method "$decoy_method" \\
        --enzyme "$decoy_enzyme" \\
        --config_file $protein_decoy_config \\
        --input_database $f \\
        --decoy_prefix "$decoy_prefix" \\
        --output_database decoy_database.fa
    """
}

