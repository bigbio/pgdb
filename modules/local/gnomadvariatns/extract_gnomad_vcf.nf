/**
 * Extract gnomAD VCF
 */
process EXTRACT_GNOMAD_VCF {

    container "nfcore/pgdb:1.0.0"
    
    when:
    params.gnomad

    input:
    file g

    output:
    path "*.vcf" ,emit: gnomad_vcf_files

    script:
    """
    zcat $g > ${g}.vcf
    """
}