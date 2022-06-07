/**
 * Create the decoy database using DecoyPYrat
 * Decoy sequences will have "DECOY_" prefix tag to the protein accession.
 */
process DECOY {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk' }"
    
    publishDir "${params.outdir}/", mode: params.publish_dir_mode,
        saveAs: { filename -> params.final_database_protein }

    when:
    params.decoy

    input:
    file f
    file protein_decoy_config

    output:
    path 'decoy_database.fa' ,emit: fasta_decoy_db_ch

    script:
    """
    pypgatk_cli.py generate-decoy \\
        --method "$params.decoy_method" \\
        --enzyme "$params.decoy_enzyme" \\
        --config_file $params.protein_decoy_config \\
        --input_database $f \\
        --decoy_prefix "$params.decoy_prefix" \\
        --output_database decoy_database.fa
    """
}

