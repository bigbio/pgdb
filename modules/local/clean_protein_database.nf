/**
 * clean the database for stop codons, and unwanted AA like: *, also remove proteins with less than 6 AA
 */
process CLEAN_PROTEIN_DATABASE {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"

    publishDir "${params.outdir}/", mode: params.publish_dir_mode,
        // Final step if not creating a decoy database - save output to params.final_database_protein
        saveAs: { filename ->
            params.decoy ? null : params.final_database_protein
        }

    when:
    params.clean_database

    input:
    file file
    file ensembl_config
    val minimum_aa

    output:
    path 'database_clean.fa' ,emit: clean_database_sh

    script:
    stop_codons = ''
    if (params.add_stop_codons){
        stop_codons = "--add_stop_codons"
    }

    """
    pypgatk_cli.py ensembl-check \\
        -in "$file" \\
        --config_file $ensembl_config \\
        -out database_clean.fa \\
        --num_aa "$minimum_aa" \\
        "$stop_codons"
    """
}
