/**
 * Download data from ensembl for the particular species.
 */
process ENSEMBL_FASTA_DOWNLOAD {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"
    
    when:
    params.add_reference ||  params.ensembl || params.altorfs || params.ncrna || params.pseudogenes || params.vcf

    input:
    file ensembl_downloader_config
    val ensembl_name

    output:
    path 'database_ensembl/*.pep.all.fa', emit: ensembl_protein_database_sub
    path 'database_ensembl/*cdna.all.fa', emit:  ensembl_cdna_database_sub
    path 'database_ensembl/*ncrna.fa', emit:  ensembl_ncrna_database_sub
    path 'database_ensembl/*.dna*.fa' ,emit: genome_fasta
    path 'database_ensembl/*.gtf' ,emit: gtf

    script:
    """
    pypgatk_cli.py ensembl-downloader \\
        --config_file $ensembl_downloader_config \\
        --ensembl_name $ensembl_name \\
        -sv -sc
    """
}
