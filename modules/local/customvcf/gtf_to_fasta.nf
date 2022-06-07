/**
 * Generate protein databse for a given VCF
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
