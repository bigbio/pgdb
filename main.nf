#!/usr/bin/env nextflow

/*
========================================================================================
                         nf-core/pgdb
========================================================================================
 nf-core/pgdb Analysensembl_downloader_configis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/pgdb
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/pgdb --ensembl_name homo_sapiens --ncrna false --pseudogenes false --altorfs false --ensembl false --gnomad false --cosmic false --cosmic_celllines false --cbioportal false --decoy false --add_reference true --vcf false

    Main arguments:
      --final_database_protein           Output file name for the final database protein fasta file under the outdir/ directory.
      --help                             Print this help document

    Process flags:
      --ncrna                            Generate protein database from non-coding RNAs [true | false] (default: false)
      --pseudogenes                      Generate protein database from pseudogenes [true | false] (default: false)
      --altorfs                          Generate alternative ORFs from canonical proteins [true | false] (default: false)
      --cbioportal                       Download cBioPortal studies and genrate protein database [true | false] (default: false)
      --cosmic                           Download COSMIC mutation files and generate protein database [true | false] (default: false)
      --cosmic_celllines                 Download COSMIC cell line files and generate protein database [true | false] (default: false)
      --ensembl                          Download ENSEMBL variants and generate protein database [true | false] (default: false)
      --gnomad                           Download gnomAD files and generate protein database [true | false] (default: false)


      --vcf                              Enable translation of a given VCF file [true | false ] (default: false)

      --add_reference                    Add the reference proteome to the file [true | false ] (default: true)

    Clean database:
      --clean_database                   Clean the database for stop codons, short protein sequences, (default: false)
      --minimum_aa                       Minimum number of AminoAcids for a protein to be included in the database (default: 6)
      --add_stop_codons                  If an stop codons is found, create two proteins from it (default: true)

    Decoy generation:
      --decoy                            Append the decoy proteins to the database [true | false] (default: false)
      --decoy_prefix                     String to be used as prefix for the generated decoy sequences
      --decoy_method                     Method used to generate the decoy database ['protein-reverse', 'protein-shuffle', 'decoypyrat'](Default: decoypyrat)
      --decoy_enzyme                     Enzyme used to generate the decoy (default: Trypsin)

    Configuration files:                 By default all config files are located in the configs directory.
      --ensembl_downloader_config        Path to configuration file for ENSEMBL download parameters
      --ensembl_config                   Path to configuration file for parameters in generating
                                         protein databases from ENSMEBL sequences
      --cosmic_config                    Path to configuration file for parameters in generating
                                         protein databases from COSMIC mutations
      --cbioportal_config                Path to configuration file for parameters in generating
                                         protein databases from cBioPortal mutations
      --protein_decoy_config             Path to configuration file for parameters used in generating
                                         decoy databases

    Database parameters:
      --taxonomy                         Taxonomy (Taxon ID) for the species to download ENSEMBL data,
                                         default is 9606 for humans. For the list of supported taxonomies see:
                                           https://www.ensembl.org/info/about/species.html

      --ensembl_name                     Ensembl Name is used to find the specific name in ENSEMBL for the taxonomy for download
                                         The list can be found here: configs/ensembl_species.txt

      --cosmic_cancer_type               Specify a tissue type to limit the COSMIC mutations to a particular caner type
                                         (by default all tumor types are used)

      --cosmic_cellline_name             Specify a sample name to limit the COSMIC cell line mutations to
                                         a particular  cell line (by default all cell lines are used)

      --cbioportal_accepted_values       Specify a tissue type to limit the cBioPortal mutations to
                                         a particular caner type (by default all tumor types are used)
      --cbioportal_filter_column         Specify a column from the clincal sample file to be used for filterring records
                                         Only values listed in cbioportal_accepted_values parameter are included, default is CANCER_TYPE
      --af_field                         Allele frequency identifier string in VCF Info column, if no AF info is given set it to empty.
                                         For human VCF files from ENSEMBL the default is set to MAF

    Data download parameters:
      --cosmic_user_name                 User name (or email) for COSMIC account
      --cosmic_password                  Password for COSMIC account
                                         In order to be able to download COSMIC data, the user should
                                         provide a user and password. Please first register in COSMIC
                                         database (https://cancer.sanger.ac.uk/cosmic/register).
      --cbioportal_study_id              Download mutations from a specific study in cbiportal
                                         default is all which downloads mutations from all studies

      --gencode_url                      URL for downloading GENCODE datafiles: gencode.v19.pc_transcripts.fa.gz and
                                         gencode.v19.annotation.gtf.gz
      --gnomad_file_url                  URL for downloading gnomAD VCF file(s)

      --vcf_file                         VCF file path to be translated
                                         Generate variants proteins by modifying sequences of affected transcripts.
                                         In case of already annotated variants it only considers variants within
                                         potential coding regions of the transcript (CDSs & stop codons for protein-coding genes, exons for non-protein coding genes)
                                         In case of not annotated variants, it considers all variants overlapping CDSs

    Output parameters:
      --publish_dir_mode [str]           Mode for publishing results in the output directory. Available:
                                         symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --outdir                           Output folder for the results by default is $baseDir/result
      --email [email]                    Set this parameter to your e-mail address to get a summary e-mail with
                                         details of the run sent to you when the workflow exits
      --email_on_fail [email]            Same as --email, except only send mail if the workflow is not successful
      -name [str]                        Name for the pipeline run. If not specified, Nextflow will automatically generate a random

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)
ensembl_downloader_config = file(params.ensembl_downloader_config, checkIfExists: true)
ensembl_config = file(params.ensembl_config)
cosmic_config = file(params.cosmic_config)
cbioportal_config = file(params.cbioportal_config)
protein_decoy_config = file(params.protein_decoy_config)

params.cbioportal_study_id = "all"

af_field = params.af_field
if (params.ensembl_name == "homo_sapiens"){
	af_field = "MAF"
}

// Pipeline checks
if ((params.cosmic || params.cosmic_celllines) && (params.cosmic_user_name=="" || params.cosmic_password=="")){
	exit 1, "User name and password has to be provided. In order to be able to download COSMIC data, the user should provide a user and password. Please first register in COSMIC database (https://cancer.sanger.ac.uk/cosmic/register)."
}

// Pipeline OS-specific commands
ZCAT = (System.properties['os.name'] == 'Mac OS X' ? 'gzcat' : 'zcat')


/**
 * Download data from ensembl for the particular species.
 */
process ensembl_fasta_download{

   when:
   params.add_reference ||  params.ensembl || params.altorfs || params.ncrna || params.pseudogenes || params.vcf

   input:
   file ensembl_downloader_config

   output:
   file "database_ensembl/*.gz" into ensembl_fasta_gz_databases

   script:
   """
   pypgatk_cli.py ensembl-downloader --config_file ${ensembl_downloader_config} --ensembl_name ${params.ensembl_name} -sv -sc
   gzip -t database_ensembl/*.gz
   """
}

/**
 * Decompress all the data downloaded from ENSEMBL
 */
process gunzip_ensembl_files{

   publishDir "${params.outdir}", mode: 'copy', overwrite: true

   input:
   file(fasta_file) from ensembl_fasta_gz_databases

   output:
   file '*.pep.all.fa' into ensembl_protein_database_sub
   file '*cdna.all.fa' into ensembl_cdna_database, ensembl_cdna_database_sub
   file '*ncrna.fa' into ensembl_ncrna_database, ensembl_ncrna_database_sub
   file '*.dna*.fa' into genome_fasta
   file '*.gtf' into gtf

   script:
   """
   gunzip -d -f ${fasta_file}
   """
}

process add_reference_proteome{

   when:
   params.add_reference

   input:
   file reference_proteome from ensembl_protein_database_sub

   output:
   file 'reference_proteome.fa' into ensembl_protein_database

   script:
   """
   cat ${reference_proteome} >> reference_proteome.fa
   """

}

/**
 * Concatenate cDNA and ncRNA databases
 **/
process merge_cdnas{

   input:
   file a from ensembl_cdna_database_sub.collect()
   file b from ensembl_ncrna_database_sub.collect()

   output:
   file 'total_cdnas.fa' into total_cdnas

   script:
   """
   cat ${a} >> total_cdnas.fa
   cat ${b} >> total_cdnas.fa
   """
}

/**
 * Creates the ncRNA protein database
 */
process add_ncrna{

   publishDir "${params.outdir}", mode: 'copy', overwrite: true

   when:
   params.ncrna

   input:
   file x from total_cdnas
   file ensembl_config

   output:
   file 'ncRNAs_proteinDB.fa' into optional_ncrna

   script:
   """
   pypgatk_cli.py dnaseq-to-proteindb --config_file "${ensembl_config}" --input_fasta ${x} --output_proteindb ncRNAs_proteinDB.fa --include_biotypes "${params.biotypes['ncRNA']}" --skip_including_all_cds --var_prefix ncRNA_
   """
}

merged_databases = ensembl_protein_database.mix(optional_ncrna)

/**
 * Creates the pseudogenes protein database
 */
process add_pseudogenes {

   publishDir "${params.outdir}", mode: 'copy', overwrite: true

   when:
   params.pseudogenes

   input:
   file x from total_cdnas
   file ensembl_config

   output:
   file 'pseudogenes_proteinDB.fa' into optional_pseudogenes

   script:
   """
   pypgatk_cli.py dnaseq-to-proteindb --config_file "${ensembl_config}" --input_fasta "${x}" --output_proteindb pseudogenes_proteinDB.fa --include_biotypes "${params.biotypes['pseudogene']}" --skip_including_all_cds --var_prefix pseudo_
   """
}

merged_databases = merged_databases.mix(optional_pseudogenes)

/**
 * Creates the altORFs protein database
 */
process add_altorfs {

   publishDir "${params.outdir}", mode: 'copy', overwrite: true

   when:
   params.altorfs

   input:
   file x from ensembl_cdna_database
   file ensembl_config

   output:
   file('altorfs_proteinDB.fa') into optional_altorfs

   script:
   """
   pypgatk_cli.py dnaseq-to-proteindb --config_file "${ensembl_config}" --input_fasta "${x}" --output_proteindb altorfs_proteinDB.fa --include_biotypes "${params.biotypes['protein_coding']}'" --skip_including_all_cds --var_prefix altorf_
   """
}

merged_databases = merged_databases.mix(optional_altorfs)

/* Mutations to proteinDB */

/**
 * Download COSMIC Mutations
 */
process cosmic_download {

	  when:
  	  params.cosmic || params.cosmic_celllines

	  input:
	  file cosmic_config

	  output:
    file "database_cosmic/All_COSMIC_Genes.fasta" into cosmic_genes
    file "database_cosmic/CosmicMutantExport.tsv" into cosmic_mutations
    file "database_cosmic/All_CellLines_Genes.fasta" into cosmic_celllines_genes
    file "database_cosmic/CosmicCLP_MutantExport.tsv" into cosmic_celllines_mutations

	  script:
	  """
	  pypgatk_cli.py cosmic-downloader --config_file "${cosmic_config}" --username ${params.cosmic_user_name} --password ${params.cosmic_password}
    gunzip -d -f database_cosmic/*.gz
	  """
}

/**
 * Generate proteindb from cosmic mutations
*/
process cosmic_proteindb{

	  publishDir "${params.outdir}", mode: 'copy', overwrite: true

	  when:
  	  params.cosmic

	  input:
	  file g from cosmic_genes
	  file m from cosmic_mutations
	  file cosmic_config

	  output:
	  file 'cosmic_proteinDB*.fa' into cosmic_proteindbs

	  script:
	  """
	  pypgatk_cli.py cosmic-to-proteindb --config_file "${cosmic_config}" --input_mutation ${m} --input_genes ${g} --filter_column 'Histology subtype 1' --accepted_values ${params.cosmic_cancer_type} --output_db cosmic_proteinDB.fa
	  """
}

merged_databases = merged_databases.mix(cosmic_proteindbs)

/**
 * Generate proteindb from cosmic cell lines mutations
*/
process cosmic_celllines_proteindb{

	  publishDir "${params.outdir}", mode: 'copy', overwrite: true

	  when:
  	  params.cosmic_celllines

	  input:
	  file g from cosmic_celllines_genes
	  file m from cosmic_celllines_mutations
	  file cosmic_config

	  output:
	  file 'cosmic_celllines_proteinDB*.fa' into cosmic_celllines_proteindbs

	  script:
	  """
	  pypgatk_cli.py cosmic-to-proteindb --config_file "${cosmic_config}" --input_mutation ${m} --input_genes ${g} --filter_column 'Sample name' --accepted_values ${params.cosmic_cellline_name} --output_db cosmic_celllines_proteinDB.fa
	  """
}

merged_databases = merged_databases.mix(cosmic_celllines_proteindbs)

/**
 * Download VCF files from ensembl for the particular species.
 */
process ensembl_vcf_download{

   when:
    params.ensembl

   input:
   file ensembl_downloader_config

   output:
   file "database_ensembl/*.vcf.gz" into ensembl_vcf_gz_files

   script:
   """
   pypgatk_cli.py ensembl-downloader --config_file ${ensembl_downloader_config} --ensembl_name ${params.ensembl_name} -sg -sp -sc -sd -sn
   """
}

/**
 * Decompress vcf files downloaded from ENSEMBL
 */
process gunzip_vcf_ensembl_files{

   label 'process_medium'
   label 'process_single_thread'

   when:
    params.ensembl

   input:
   file vcf_file from ensembl_vcf_gz_files.flatten().map{ file(it) }

   output:
   file "*.vcf" into ensembl_vcf_files

   script:
   """
   gunzip -d -f $vcf_file
   """
}

process check_ensembl_vcf{

   label 'process_medium'
   label 'process_single_thread'

   when:
   params.ensembl

   input:
   file vcf_file from ensembl_vcf_files

   output:
   file "checked_*.vcf" into ensembl_vcf_files_checked

   script:
   """
   awk 'BEGIN{FS=OFS="\t"}{if(\$1~"#" || (\$5!="" && \$4!="")) print}' $vcf_file > checked_$vcf_file
   """
}

/**
 * Generate protein database(s) from ENSEMBL vcf file(s)
 */
process ensembl_vcf_proteinDB {

   label 'process_medium'
   label 'process_single_thread'

   when:
   params.ensembl

   input:
   file v from ensembl_vcf_files_checked
   file f from total_cdnas
   file g from gtf
   file e from ensembl_config

   output:
   file "${v}_proteinDB.fa" into proteinDB_vcf

   script:
   """
   pypgatk_cli.py vcf-to-proteindb --config_file ${e} --af_field "${af_field}" --input_fasta ${f} --gene_annotations_gtf ${g} --vcf ${v} --output_proteindb "${v}_proteinDB.fa"  --var_prefix ensvar --annotation_field_name 'CSQ'
   """
}

//concatenate all ensembl proteindbs into one
proteinDB_vcf
	.collectFile(name: 'ensembl_proteindb.fa', newLine: false, storeDir: "${baseDir}/result")
	.set {proteinDB_vcf_final}

merged_databases = merged_databases.mix(proteinDB_vcf_final)

/****** Custom VCF      *****/
/**
 * Generate protein databse for a given VCF
 */
process gtf_to_fasta {

   when:
   params.vcf

   input:
   file g from gtf
   file f from genome_fasta

   output:
   file "transcripts.fa" into gtf_transcripts_fasta

   script:
   """
   gffread -w transcripts.fa -g ${f} ${g}
   """
}
process vcf_proteinDB {

   publishDir "${params.outdir}", mode: 'copy', overwrite: true

   when:
   params.vcf

   input:
   file v from params.vcf_file
   file f from gtf_transcripts_fasta
   file g from gtf
   file e from ensembl_config

   output:
   file "${v}_proteinDB.fa" into proteinDB_custom_vcf

   script:
   """
   pypgatk_cli.py vcf-to-proteindb --config_file ${e} --af_field "${af_field}" --input_fasta ${f} --gene_annotations_gtf ${g} --vcf ${params.vcf_file} --output_proteindb ${params.vcf_file.baseName}_proteinDB.fa --annotation_field_name ''
   """
}

merged_databases = merged_databases.mix(proteinDB_custom_vcf)


/****** gnomAD variatns *****/

/**
 * Download gencode files (fasta and gtf)
 */
process gencode_download{

   when:
	  params.gnomad

   input:
	 val g from params.gencode_url

   output:
	 file("gencode.v19.pc_transcripts.fa") into gencode_fasta
	 file("gencode.v19.annotation.gtf") into gencode_gtf

   script:
	 """
	 wget ${g}/gencode.v19.pc_transcripts.fa.gz
	 wget ${g}/gencode.v19.annotation.gtf.gz
	 gunzip *.gz
	 """
}

/**
 * Download gnomAD variants (VCF) - requires gsutil
 */
process gnomad_download{

   when:
	  params.gnomad

   input:
	 val g from params.gnomad_file_url

   output:
   file "*.vcf.bgz" into gnomad_vcf_bgz

   script:
   """
   gsutil cp ${g} .
   """
}

/**
 * Extract gnomAD VCF
 */
process extract_gnomad_vcf{

   when:
   params.gnomad

   input:
   file g from gnomad_vcf_bgz.flatten().map{ file(it) }

   output:
   file "*.vcf" into gnomad_vcf_files

   script:
   """
   zcat ${g} > ${g}.vcf
   """
}

/**
 * Generate gmomAD proteinDB
 */
process gnomad_proteindb{

   when:
   params.gnomad

   input:
   file v from gnomad_vcf_files
   file f from gencode_fasta
   file g from gencode_gtf
   file e from ensembl_config

   output:
   file "${v}_proteinDB.fa" into gnomad_vcf_proteindb

   script:
   """
   pypgatk_cli.py vcf-to-proteindb --config_file ${e} --vcf ${v} --input_fasta ${f} --gene_annotations_gtf ${g} --output_proteindb "${v}_proteinDB.fa" --af_field controls_AF --transcript_index 6 --annotation_field_name vep  --var_prefix gnomadvar
   """
}

//concatenate all gnomad proteindbs into one
gnomad_vcf_proteindb
	.collectFile(name: 'gnomad_proteindb.fa', newLine: false, storeDir: "${baseDir}/result")
	.set {gnomad_vcf_proteindb_final}

merged_databases = merged_databases.mix(gnomad_vcf_proteindb_final)

/****** cBioPortal mutations *****/
/**
 * Download GRCh37 CDS file from ENSEMBL release 75
 */
process cds_GRCh37_download{

   when:
   params.cbioportal

   output:
   file("Homo_sapiens.GRCh37.75.cds.all.fa") into GRCh37_cds

   script:
   """
   wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.75.cds.all.fa.gz
   gunzip *.gz
   """
}

/**
 * Download all cBioPortal studies using git-lfs
*/
 process download_all_cbioportal {

    when:
         params.cbioportal

   output:
 	 file('cbioportal_allstudies_data_mutations_mskcc.txt') into cbio_mutations
 	 file('cbioportal_allstudies_data_clinical_sample.txt') into cbio_samples

   script:
   if (params.cbioportal_study_id == "all")
        """
        git clone https://github.com/cBioPortal/datahub.git
        cd datahub
        git lfs install --local --skip-smudge
        git lfs pull -I public --include "data*clinical*sample.txt"
        git lfs pull -I public --include "data_mutations_mskcc.txt"
        cd ..
        cat datahub/public/*/data_mutations_mskcc.txt > cbioportal_allstudies_data_mutations_mskcc.txt
        cat datahub/public/*/*data*clinical*sample.txt | awk 'BEGIN{FS=OFS="\\t"}{if(\$1!~"#SAMPLE_ID"){gsub("#SAMPLE_ID", "\\nSAMPLE_ID");} print}' | awk 'BEGIN{FS=OFS="\\t"}{s=0; j=0; for(i=1;i<=NF;i++){if(\$i=="CANCER_TYPE_DETAILED") j=1; if(\$i=="CANCER_TYPE") s=1;} if(j==1 && s==0){gsub("CANCER_TYPE_DETAILED", "CANCER_TYPE");} print;}' > cbioportal_allstudies_data_clinical_sample.txt
        """
    else
        """
        pypgatk_cli.py cbioportal-downloader --config_file "${cbioportal_config}" -d "${params.cbioportal_study_id}"
        tar -xzvf database_cbioportal/${params.cbioportal_study_id}.tar.gz
        cat ${params.cbioportal_study_id}/data_mutations_mskcc.txt > cbioportal_allstudies_data_mutations_mskcc.txt
        cat ${params.cbioportal_study_id}/data_clinical_sample.txt | awk 'BEGIN{FS=OFS="\\t"}{if(\$1!~"#SAMPLE_ID"){gsub("#SAMPLE_ID", "\\nSAMPLE_ID");} print}' | awk 'BEGIN{FS=OFS="\\t"}{s=0; j=0; for(i=1;i<=NF;i++){if(\$i=="CANCER_TYPE_DETAILED") j=1; if(\$i=="CANCER_TYPE") s=1;} if(j==1 && s==0){gsub("CANCER_TYPE_DETAILED", "CANCER_TYPE");} print;}' > cbioportal_allstudies_data_clinical_sample.txt
        """
 }

/**
 * Generate proteinDB from cBioPortal mutations
 */
 process cbioportal_proteindb{

   publishDir "${params.outdir}", mode: 'copy', overwrite: true

   when:
   params.cbioportal

   input:
   file g from GRCh37_cds
   file m from cbio_mutations
   file s from cbio_samples
   file cbioportal_config

   output:
   file 'cbioPortal_proteinDB*.fa' into cBioportal_proteindb

   script:
   """
   pypgatk_cli.py cbioportal-to-proteindb --config_file ${cbioportal_config} --input_mutation ${m} --input_cds ${g} --clinical_sample_file ${s} --filter_column ${params.cbioportal_filter_column} --accepted_values ${params.cbioportal_accepted_values} --output_db cbioPortal_proteinDB.fa
   """
}

merged_databases = merged_databases.mix(cBioportal_proteindb)

/**
 * Concatenate all generated databases from merged_databases channel to the final_database_protein file
 */
process merge_proteindbs {

   publishDir "${params.outdir}", mode: 'copy', overwrite: true

   input:
   file("proteindb*") from merged_databases.collect()

   output:
   file 'merged_databases.fa' into to_clean_ch

   script:
   """
   cat proteindb* > merged_databases.fa
   """
}

stop_codons = ''
if (params.add_stop_codons){
	stop_codons = "--add_stop_codons"
}

/**
 * clean the database for stop codons, and unwanted AA like: *, also remove proteins with less than 6 AA
 */
process clean_protein_database {

   publishDir "${params.outdir}", mode: 'copy', overwrite: true

   when:
     params.clean_database

   input:
   file file from to_clean_ch
   file e from ensembl_config

   output:
   file 'database_clean.fa' into clean_database_sh

   script:
   """
   pypgatk_cli.py ensembl-check -in "${file}" --config_file "${e}" -out database_clean.fa --num_aa "${params.minimum_aa}" "${stop_codons}"
   """
}

to_protein_decoy_ch = params.clean_database ? clean_database_sh : to_clean_ch

/**
 * Create the decoy database using DecoyPYrat
 * Decoy sequences will have "DECOY_" prefix tag to the protein accession.
 */
process decoy {

   publishDir "${params.outdir}", mode: 'copy', overwrite: true

   when:
    params.decoy

   input:
   file f from to_protein_decoy_ch
   file protein_decoy_config

   output:
   file 'decoy_database.fa' into fasta_decoy_db_ch

   script:
   """
   pypgatk_cli.py generate-decoy --method "${params.decoy_method}" --enzyme "${params.decoy_enzyme}" --config_file ${protein_decoy_config} --input_database $f --decoy_prefix "${params.decoy_prefix}" --output_database decoy_database.fa
   """
}

result_database_ch = params.decoy ? fasta_decoy_db_ch: to_protein_decoy_ch

/** Write the final results to S3 bucket**/

result_database_ch.subscribe { results -> results.copyTo("${params.outdir}/${params.final_database_protein}")}


//--------------------------------------------------------------- //
//---------------------- Nextflow specifics --------------------- //
//--------------------------------------------------------------- //


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
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
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-proteomicslfq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/proteomicslfq Workflow Summary'
    section_href: 'https://github.com/nf-core/proteomicslfq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {

    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/pgdb] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/pgdb] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
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

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report
    //def mqc_report = null
    //try {
    //    if (workflow.success) {
    //        mqc_report = ch_multiqc_report.getVal()
    //        if (mqc_report.getClass() == ArrayList) {
    //            log.warn "[nf-core/pgdb] Found multiple reports from process 'multiqc', will use only one"
    //            mqc_report = mqc_report[0]
    //        }
    //    }
    //} catch (all) {
    //    log.warn "[nf-core/pgdb] Could not attach MultiQC report to summary email"
    //}

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


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/pgdb v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
