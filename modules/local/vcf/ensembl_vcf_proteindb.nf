/**
 * Generate protein database(s) from ENSEMBL vcf file(s)
 */
process ENSEMBL_VCF_PROTEINDB {

    label 'process_medium'
    label 'process_single_thread'

    container "nfcore/pgdb:1.0.0"
    
    when:
    params.ensembl

    input:
    file v
    file f
    file g
    file e
    val ensembl_af_field 

    output:
    path "${v}_proteinDB.fa" ,emit: proteinDB_vcf

    script:
    """
    pypgatk_cli.py vcf-to-proteindb \\
        --config_file $e \\
        --af_field $ensembl_af_field \\
        --input_fasta $f \\
        --gene_annotations_gtf $g \\
        --vcf $v \\
        --output_proteindb "${v}_proteinDB.fa"  \\
        --var_prefix ensvar \\
        --annotation_field_name 'CSQ'
    """
}