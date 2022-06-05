/**
 * Download data from ensembl for the particular species.
 */
process ENSEMBL_FASTA_DOWNLOAD {

    container "nfcore/pgdb:1.0.0"
    
    when:
    params.add_reference ||  params.ensembl || params.altorfs || params.ncrna || params.pseudogenes || params.vcf

    input:
    file ensembl_downloader_config

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
        --ensembl_name $params.ensembl_name \\
        -sv -sc
    """
}