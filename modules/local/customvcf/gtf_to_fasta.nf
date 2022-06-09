/**
 * Generate protein databse for a given VCF
 conda (params.enable_conda ? "bioconda::gffread=0.12.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread_0.12.1--h8b12597_0' :
        'quay.io/biocontainers/pypgatk:0.12.1--h8b12597_0' }"
 */
process GTF_TO_FASTA {

    
    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"
    
    when:
    params.vcf

    input:
    file g
    file f

    output:
    path "transcripts.fa" ,emit:gtf_transcripts_fasta

    script:
    """
    gffread -w transcripts.fa -g $f $g
    """
}
