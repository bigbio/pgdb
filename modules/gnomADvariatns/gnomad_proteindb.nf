/**
 * Generate gmomAD proteinDB
 */
process GNOMAD_PROTEINDB {

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