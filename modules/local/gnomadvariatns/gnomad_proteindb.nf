/**
 * Generate gmomAD proteinDB
 */
process GNOMAD_PROTEINDB {

    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"
    
    when:
    params.gnomad

    input:
    file v
    file f
    file g
    file e

    output:
    path "${v}_proteinDB.fa" ,emit:gnomad_vcf_proteindb

    script:
    """
    pypgatk_cli.py vcf-to-proteindb \\
        --config_file $e \\
        --vcf $v \\
        --input_fasta $f \\
        --gene_annotations_gtf $g \\
        --output_proteindb "${v}_proteinDB.fa" \\
        --af_field controls_AF \\
        --transcript_index 6 \\
        --annotation_field_name vep  \\
        --var_prefix gnomadvar
    """
}
