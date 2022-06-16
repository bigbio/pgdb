#!/usr/bin/env nextflow

/*
========================================================================================
                         nf-core/pgdb
========================================================================================
 nf-core/pgdb Proteogenomics database generation
 #### Homepage / Documentation
 https://github.com/nf-core/pgdb
----------------------------------------------------------------------------------------
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPgdb.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// if (params.validate_params) {
//     NfcoreSchema.validateParameters(params, log)
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
ensembl_downloader_config = file(params.ensembl_downloader_config, checkIfExists: true)
ensembl_config = file(params.ensembl_config)
cosmic_config = file(params.cosmic_config)
if (params.cosmicgenes&&params.cosmicmutations) {
    cosmicgenes = file(params.cosmicgenes)
    cosmicmutations = file(params.cosmicmutations)
}
if (params.cosmiccelllines_genes&&params.cosmiccelllines_mutations) {
    cosmiccelllines_genes = file(params.cosmiccelllines_genes)
    cosmiccelllines_mutations = file(params.cosmiccelllines_mutations)
}
cbioportal_config = file(params.cbioportal_config)
protein_decoy_config = file(params.protein_decoy_config)

af_field = params.af_field
if (af_field == null ) {
    af_field = ""
}
if (params.ensembl_name == "homo_sapiens"){
    af_field = "MAF"
}
// ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
// ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
//def modules = params.modules.clone()

include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'

include { ENSEMBL_FASTA_DOWNLOAD } from '../modules/local/proteomes/ensembl_fasta_download'
include { ADD_REFERENCE_PROTEOME } from '../modules/local/proteomes/add_reference_proteome'
include { MERGE_CDNAS} from '../modules/local/proteomes/merge_cdnas'
include { ADD_NCRNA } from '../modules/local/proteomes/add_ncrna'
include { ADD_PSEUDOGENES } from '../modules/local/proteomes/add_pseudogenes'
include { ADD_ALTORFS } from '../modules/local/proteomes/add_altorfs'

include { COSMIC_DOWNLOAD } from '../modules/local/cosmicmutations/cosmic_download'
include { COSMIC_PROTEINDB } from '../modules/local/cosmicmutations/cosmic_proteindb'
include { COSMIC_CELLLINES_PROTEINDB } from '../modules/local/cosmicmutations/cosmic_celllines_proteindb'
include { COSMIC_PROTEINDB_LOCAL } from '../modules/local/cosmicmutations/cosmic_proteindb_local'
include { COSMIC_CELLLINES_PROTEINDB_LOCAL } from '../modules/local/cosmicmutations/cosmic_celllines_proteindb_local'

include { ENSEMBL_VCF_DOWNLOAD } from '../modules/local/vcf/ensembl_vcf_download'
include { CHECK_ENSEMBL_VCF } from '../modules/local/vcf/check_ensembl_vcf'
include { ENSEMBL_VCF_PROTEINDB } from '../modules/local/vcf/ensembl_vcf_proteindb'

include { GTF_TO_FASTA } from '../modules/local/customvcf/gtf_to_fasta'
include { VCF_PROTEINDB } from '../modules/local/customvcf/vcf_proteinDB'

include { GENCODE_DOWNLOAD } from '../modules/local/gnomadvariatns/gencode_download'
include { GNOMAD_DOWNLOAD } from '../modules/local/gnomadvariatns/gnomad_download'
include { EXTRACT_GNOMAD_VCF } from '../modules/local/gnomadvariatns/extract_gnomad_vcf'
include { GNOMAD_PROTEINDB } from '../modules/local/gnomadvariatns/gnomad_proteindb'

include { CDS_GRCH37_DOWNLOAD } from '../modules/local/cbioportalmutations/cds_GRCh37_download'
include { DOWNLOAD_ALL_CBIOPORTAL } from '../modules/local/cbioportalmutations/download_all_cbioportal'
include { CBIOPORTAL_PROTEINDB } from '../modules/local/cbioportalmutations/cbioportal_proteindb'

include { MERGE_PROTEINDBS } from '../modules/local/merge_proteindbs'
include { CLEAN_PROTEIN_DATABASE } from '../modules/local/clean_protein_database'
include { DECOY } from '../modules/local/decoy'
include { OUTPUT_DOCUMENTATION } from '../modules/local/output_documentation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PGDB {

//     // Parse software version numbers
//     GET_SOFTWARE_VERSIONS()

    // Download data from ensembl for the particular species
    ENSEMBL_FASTA_DOWNLOAD(ensembl_downloader_config,params.ensembl_name)

    ADD_REFERENCE_PROTEOME(ENSEMBL_FASTA_DOWNLOAD.out.ensembl_protein_database_sub)

    //Concatenate cDNA and ncRNA databases
    MERGE_CDNAS(ENSEMBL_FASTA_DOWNLOAD.out.ensembl_cdna_database_sub.collect(),ENSEMBL_FASTA_DOWNLOAD.out.ensembl_ncrna_database_sub.collect())

    //Creates the ncRNA protein database
    ADD_NCRNA(MERGE_CDNAS.out.total_cdnas,ensembl_config)
    merged_databases = ADD_REFERENCE_PROTEOME.out.ensembl_protein_database.mix(ADD_NCRNA.out.optional_ncrna)

    //Creates the pseudogenes protein database
    ADD_PSEUDOGENES(MERGE_CDNAS.out.total_cdnas,ensembl_config)
    merged_databases = merged_databases.mix(ADD_PSEUDOGENES.out.optional_pseudogenes)

    //Creates the altORFs protein database
    ADD_ALTORFS(ENSEMBL_FASTA_DOWNLOAD.out.ensembl_cdna_database_sub,ensembl_config)
    merged_databases = merged_databases.mix(ADD_ALTORFS.out.optional_altorfs)


    /* Mutations to proteinDB */

    //Download COSMIC Mutations
    COSMIC_DOWNLOAD()

    //Generate proteindb from cosmic mutations
    COSMIC_PROTEINDB(COSMIC_DOWNLOAD.out.cosmic_genes,COSMIC_DOWNLOAD.out.cosmic_mutations,cosmic_config,params.cosmic_cancer_type)
    if (params.cosmic) {
        merged_databases = merged_databases.mix(COSMIC_PROTEINDB.out.cosmic_proteindbs)
    }

    //Generate proteindb from local cosmic mutations
    if (params.cosmicgenes&&params.cosmicmutations) {
        COSMIC_PROTEINDB_LOCAL(cosmic_config,params.cosmic_cancer_type,cosmicmutations,cosmicgenes)
        merged_databases = merged_databases.mix(COSMIC_PROTEINDB_LOCAL.out.cosmic_proteindbs_uselocal)
    }

    //Generate proteindb from cosmic cell lines mutations
    COSMIC_CELLLINES_PROTEINDB(COSMIC_DOWNLOAD.out.cosmic_celllines_genes,COSMIC_DOWNLOAD.out.cosmic_celllines_mutations,cosmic_config,params.cosmic_cellline_name)
    if (params.cosmic_celllines) {
        merged_databases = merged_databases.mix(COSMIC_CELLLINES_PROTEINDB.out.cosmic_celllines_proteindbs)
    }

    //Generate proteindb from local cosmic cell lines mutations
    if (params.cosmiccelllines_genes&&params.cosmiccelllines_mutations) {
        COSMIC_CELLLINES_PROTEINDB_LOCAL(cosmic_config,params.cosmic_cellline_name,cosmiccelllines_mutations,cosmiccelllines_genes)
        merged_databases = merged_databases.mix(COSMIC_CELLLINES_PROTEINDB_LOCAL.out.cosmic_celllines_proteindbs_uselocal)
    }

    //Download VCF files from ensembl for the particular species
    ENSEMBL_VCF_DOWNLOAD(ensembl_downloader_config,params.ensembl_name)
    CHECK_ENSEMBL_VCF(ENSEMBL_VCF_DOWNLOAD.out.ensembl_vcf_files)

    //Generate protein database(s) from ENSEMBL vcf file(s)
    ENSEMBL_VCF_PROTEINDB(CHECK_ENSEMBL_VCF.out.ensembl_vcf_files_checked,MERGE_CDNAS.out.total_cdnas,ENSEMBL_FASTA_DOWNLOAD.out.gtf,ensembl_config,af_field)

    //concatenate all ensembl proteindbs into one
    ENSEMBL_VCF_PROTEINDB.out.proteinDB_vcf.collectFile(name: 'ensembl_proteindb.fa', newLine: false, storeDir: "${projectDir}/result")
        .set {proteinDB_vcf_final}
    merged_databases = merged_databases.mix(proteinDB_vcf_final)

    /*Custom VCF */

    //Generate protein databse for a given VCF
    GTF_TO_FASTA(ENSEMBL_FASTA_DOWNLOAD.out.gtf,ENSEMBL_FASTA_DOWNLOAD.out.genome_fasta)
    vcf_file = params.vcf_file ? Channel.fromPath(params.vcf_file, checkIfExists: true) : Channel.empty()

    VCF_PROTEINDB(vcf_file,GTF_TO_FASTA.out.gtf_transcripts_fasta,ENSEMBL_FASTA_DOWNLOAD.out.gtf,ensembl_config,af_field)
    merged_databases = merged_databases.mix(VCF_PROTEINDB.out.proteinDB_custom_vcf)

    /*gnomAD variatns */

    //Download gencode files (fasta and gtf)
    GENCODE_DOWNLOAD(params.gencode_url)

    //Download gnomAD variants (VCF) - requires gsutil
    GNOMAD_DOWNLOAD(params.gnomad_file_url)

    //Extract gnomAD VCF
    EXTRACT_GNOMAD_VCF(GNOMAD_DOWNLOAD.out.gnomad_vcf_bgz.flatten().map{ file(it) })

    //Generate gmomAD proteinDB
    GNOMAD_PROTEINDB(EXTRACT_GNOMAD_VCF.out.gnomad_vcf_files,GENCODE_DOWNLOAD.out.gencode_fasta,GENCODE_DOWNLOAD.out.gencode_gtf,params.ensembl_config)

    //concatenate all gnomad proteindbs into one
    GNOMAD_PROTEINDB.out.gnomad_vcf_proteindb.collectFile(name: 'gnomad_proteindb.fa', newLine: false, storeDir: "${projectDir}/result")
        .set {gnomad_vcf_proteindb_final}
    merged_databases = merged_databases.mix(gnomad_vcf_proteindb_final)

    /*cBioPortal mutations */

    //Download GRCh37 CDS file from ENSEMBL release 75
    CDS_GRCH37_DOWNLOAD()

    //Download all cBioPortal studies using git-lfs
    DOWNLOAD_ALL_CBIOPORTAL(cbioportal_config,params.cbioportal_study_id)

    //Generate proteinDB from cBioPortal mutations
    CBIOPORTAL_PROTEINDB(CDS_GRCH37_DOWNLOAD.out.ch_GRCh37_cds,DOWNLOAD_ALL_CBIOPORTAL.out.cbio_mutations,DOWNLOAD_ALL_CBIOPORTAL.out.cbio_samples,cbioportal_config,params.cbioportal_filter_column,params.cbioportal_accepted_values)
    merged_databases = merged_databases.mix(CBIOPORTAL_PROTEINDB.out.cBioportal_proteindb)

    //Concatenate all generated databases from merged_databases channel to the final_database_protein file
    MERGE_PROTEINDBS(merged_databases.collect())

    //clean the database for stop codons, and unwanted AA like: *, also remove proteins with less than 6 AA
    CLEAN_PROTEIN_DATABASE(MERGE_PROTEINDBS.out.to_clean_ch,ensembl_config,params.minimum_aa)

    to_protein_decoy_ch = params.clean_database ? CLEAN_PROTEIN_DATABASE.out.clean_database_sh : MERGE_PROTEINDBS.out.to_clean_ch

    //Create the decoy database using DecoyPYrat
    //Decoy sequences will have "DECOY_" prefix tag to the protein accession
    DECOY(to_protein_decoy_ch,protein_decoy_config,params.decoy_method,params.decoy_enzyme,params.decoy_prefix)

//     //Output Description HTML
//     OUTPUT_DOCUMENTATION(ch_output_docs,ch_output_docs_images)

}

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
