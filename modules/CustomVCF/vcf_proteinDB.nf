process VCF_PROTEINDB {

    when:
    params.vcf

    input:
    file v
    file f
    file g
    file e
    val af_field

    output:
    path "*_proteinDB.fa" ,emit: proteinDB_custom_vcf

    script:
    """
    awk 'BEGIN{FS=OFS="\t"}{if(\$1=="chrM") \$1="MT"; gsub("chr","",\$1); print}' \\
        $v > ${v.baseName}_changedChrNames.vcf

    pypgatk_cli.py vcf-to-proteindb \\
        --config_file $e \\
        --af_field "$af_field" \\
        --input_fasta $f \\
        --gene_annotations_gtf $g \\
        --vcf ${v.baseName}_changedChrNames.vcf \\
        --output_proteindb ${v.baseName}_proteinDB.fa \\
        --annotation_field_name ''
    """
}