/**
 * Check VCF files from ensembl for the particular species
 */
process CHECK_ENSEMBL_VCF {

    label 'process_medium'
    label 'process_single_thread'
    
    container "nfcore/pgdb:1.0.0"

    when:
    params.ensembl

    input:
    file vcf_file

    output:
    path "checked_*.vcf" ,emit:ensembl_vcf_files_checked

    script:
    """
    awk 'BEGIN{FS=OFS="\t"}{if(\$1~"#" || (\$5!="" && \$4!="")) print}' $vcf_file > checked_$vcf_file
    """
}