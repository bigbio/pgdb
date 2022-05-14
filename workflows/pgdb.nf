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

log.info Headers.nf_core(workflow, params.monochrome_logs)

/*
========================================================================================
    PRINT HELP
========================================================================================
*/
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/pgdb -profile docker --ensembl_name homo_sapiens"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

/*
========================================================================================
    VALIDATE PARAMETERS
========================================================================================
*/

if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)
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
ensembl_af_field = params.af_field
if (params.ensembl_name == "homo_sapiens"){
    ensembl_af_field = "MAF"
}

// Pipeline checks
if ((params.cosmic || params.cosmic_celllines) && (!params.cosmic_user_name || !params.cosmic_password)){
    exit 1, "User name and password has to be provided. In order to be able to download COSMIC data. Please first register in COSMIC database (https://cancer.sanger.ac.uk/cosmic/register)."
}
if ((params.cosmic&&params.cosmicgenes&&params.cosmicmutations)||(params.cosmic_celllines&&params.cosmiccelllines_genes&&params.cosmiccelllines_mutations)) {
    exit 1, "You can only choose to download data or use local data."
}
if ((params.cosmicgenes&&!params.cosmicmutations) || (!params.cosmicgenes&&params.cosmicmutations)){
    exit 1, "You have to provide both genes and mutations."
}
if ((params.cosmiccelllines_genes&&!params.cosmiccelllines_mutations) || (!params.cosmiccelllines_genes&&params.cosmiccelllines_mutations)){
    exit 1, "You have to provide both genes and mutations."
}

/*
========================================================================================
    PRINT PARAMETER SUMMARY 
========================================================================================
*/
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-pgdb-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/pgdb Workflow Summary'
    section_href: 'https://github.com/nf-core/pgdb'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }


/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
//def modules = params.modules.clone()

include { GET_SOFTWARE_VERSIONS } from '../modules/get_software_versions'

include { ENSEMBL_FASTA_DOWNLOAD } from '../modules/CanonicalAndNoncanonicalProteomes/ensembl_fasta_download'
include { ADD_REFERENCE_PROTEOME } from '../modules/CanonicalAndNoncanonicalProteomes/add_reference_proteome'
include { MERGE_CDNAS} from '../modules/CanonicalAndNoncanonicalProteomes/merge_cdnas'
include { ADD_NCRNA } from '../modules/CanonicalAndNoncanonicalProteomes/add_ncrna'
include { ADD_PSEUDOGENES } from '../modules/CanonicalAndNoncanonicalProteomes/add_pseudogenes'
include { ADD_ALTORFS } from '../modules/CanonicalAndNoncanonicalProteomes/add_altorfs'

include { COSMIC_DOWNLOAD } from '../modules/CosmicMutations/cosmic_download'
include { COSMIC_PROTEINDB } from '../modules/CosmicMutations/cosmic_proteindb'
include { COSMIC_CELLLINES_PROTEINDB } from '../modules/CosmicMutations/cosmic_celllines_proteindb'
include { COSMIC_PROTEINDB_LOCAL } from '../modules/CosmicMutations/cosmic_proteindb_local'
include { COSMIC_CELLLINES_PROTEINDB_LOCAL } from '../modules/CosmicMutations/cosmic_celllines_proteindb_local'

include { ENSEMBL_VCF_DOWNLOAD } from '../modules/VCF/ensembl_vcf_download'
include { CHECK_ENSEMBL_VCF } from '../modules/VCF/check_ensembl_vcf'
include { ENSEMBL_VCF_PROTEINDB } from '../modules/VCF/ensembl_vcf_proteindb'

include { GTF_TO_FASTA } from '../modules/CustomVCF/gtf_to_fasta'
include { VCF_PROTEINDB } from '../modules/CustomVCF/vcf_proteinDB'

include { GENCODE_DOWNLOAD } from '../modules/gnomADvariatns/gencode_download'
include { GNOMAD_DOWNLOAD } from '../modules/gnomADvariatns/gnomad_download'
include { EXTRACT_GNOMAD_VCF } from '../modules/gnomADvariatns/extract_gnomad_vcf'
include { GNOMAD_PROTEINDB } from '../modules/gnomADvariatns/gnomad_proteindb'

include { CDS_GRCH37_DOWNLOAD } from '../modules/cBioPortalMutations/cds_GRCh37_download'
include { DOWNLOAD_ALL_CBIOPORTAL } from '../modules/cBioPortalMutations/download_all_cbioportal'
include { CBIOPORTAL_PROTEINDB } from '../modules/cBioPortalMutations/cbioportal_proteindb'

include { MERGE_PROTEINDBS } from '../modules/merge_proteindbs'
include { CLEAN_PROTEIN_DATABASE } from '../modules/clean_protein_database'
include { DECOY } from '../modules/decoy'
include { OUTPUT_DOCUMENTATION } from '../modules/output_documentation'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow PGDB {

    // Parse software version numbers
    GET_SOFTWARE_VERSIONS()

    //Download data from ensembl for the particular species
    ENSEMBL_FASTA_DOWNLOAD(ensembl_downloader_config)

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
    COSMIC_DOWNLOAD(cosmic_config)

    //Generate proteindb from cosmic mutations
    COSMIC_PROTEINDB(COSMIC_DOWNLOAD.out.cosmic_genes,COSMIC_DOWNLOAD.out.cosmic_mutations,cosmic_config)
    if (params.cosmic) {
        merged_databases = merged_databases.mix(COSMIC_PROTEINDB.out.cosmic_proteindbs)
    }

    //Generate proteindb from local cosmic mutations
    if (params.cosmicgenes&&params.cosmicmutations) {
        COSMIC_PROTEINDB_LOCAL(cosmicgenes,cosmicmutations,cosmic_config)
        merged_databases = merged_databases.mix(COSMIC_PROTEINDB_LOCAL.out.cosmic_proteindbs_uselocal)
    }
    
    //Generate proteindb from cosmic cell lines mutations
    COSMIC_CELLLINES_PROTEINDB(COSMIC_DOWNLOAD.out.cosmic_celllines_genes,COSMIC_DOWNLOAD.out.cosmic_celllines_mutations,cosmic_config)
    if (params.cosmic_celllines) {
        merged_databases = merged_databases.mix(COSMIC_CELLLINES_PROTEINDB.out.cosmic_celllines_proteindbs)
    }

    //Generate proteindb from local cosmic cell lines mutations
    if (params.cosmiccelllines_genes&&params.cosmiccelllines_mutations) {
        COSMIC_CELLLINES_PROTEINDB_LOCAL(cosmiccelllines_genes,cosmiccelllines_mutations,cosmic_config)
        merged_databases = merged_databases.mix(COSMIC_CELLLINES_PROTEINDB_LOCAL.out.cosmic_celllines_proteindbs_uselocal)
    }

    //Download VCF files from ensembl for the particular species
    ENSEMBL_VCF_DOWNLOAD(ensembl_downloader_config)
    CHECK_ENSEMBL_VCF(ENSEMBL_VCF_DOWNLOAD.out.ensembl_vcf_files)

    //Generate protein database(s) from ENSEMBL vcf file(s)
    ENSEMBL_VCF_PROTEINDB(CHECK_ENSEMBL_VCF.out.ensembl_vcf_files_checked,MERGE_CDNAS.out.total_cdnas,ENSEMBL_FASTA_DOWNLOAD.out.gtf,ensembl_config,ensembl_af_field)

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
    GNOMAD_PROTEINDB(EXTRACT_GNOMAD_VCF.out.gnomad_vcf_files,GENCODE_DOWNLOAD.out.gencode_fasta,GENCODE_DOWNLOAD.out.gencode_gtf,ensembl_config)

    //concatenate all gnomad proteindbs into one
    GNOMAD_PROTEINDB.out.gnomad_vcf_proteindb.collectFile(name: 'gnomad_proteindb.fa', newLine: false, storeDir: "${projectDir}/result")
        .set {gnomad_vcf_proteindb_final}
    merged_databases = merged_databases.mix(gnomad_vcf_proteindb_final)

    /*cBioPortal mutations */

    //Download GRCh37 CDS file from ENSEMBL release 75
    CDS_GRCH37_DOWNLOAD()

    //Download all cBioPortal studies using git-lfs
    DOWNLOAD_ALL_CBIOPORTAL()

    //Generate proteinDB from cBioPortal mutations
    CBIOPORTAL_PROTEINDB(CDS_GRCH37_DOWNLOAD.out.ch_GRCh37_cds,DOWNLOAD_ALL_CBIOPORTAL.out.cbio_mutations,DOWNLOAD_ALL_CBIOPORTAL.out.cbio_samples,cbioportal_config)
    merged_databases = merged_databases.mix(CBIOPORTAL_PROTEINDB.out.cBioportal_proteindb)

    //Concatenate all generated databases from merged_databases channel to the final_database_protein file
    MERGE_PROTEINDBS(merged_databases.collect())

    //clean the database for stop codons, and unwanted AA like: *, also remove proteins with less than 6 AA
    CLEAN_PROTEIN_DATABASE(MERGE_PROTEINDBS.out.to_clean_ch,ensembl_config)

    to_protein_decoy_ch = params.clean_database ? CLEAN_PROTEIN_DATABASE.out.clean_database_sh : MERGE_PROTEINDBS.out.to_clean_ch

    //Create the decoy database using DecoyPYrat
    //Decoy sequences will have "DECOY_" prefix tag to the protein accession
    DECOY(to_protein_decoy_ch,protein_decoy_config)

    //Output Description HTML
    OUTPUT_DOCUMENTATION(ch_output_docs,ch_output_docs_images)
 
}


// Completion e-mail notification
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/pgdb] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/pgdb] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/pgdb] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
                mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/pgdb] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/pgdb]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/pgdb]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error '====================================================\n' +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            '============================================================'
                }
            }
        }
    }
}
