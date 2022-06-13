/**
 * Download all cBioPortal studies using git-lfs
*/
process DOWNLOAD_ALL_CBIOPORTAL {
    
    conda (params.enable_conda ? "bioconda::pypgatk=0.0.19 conda-forge::git-lfs=2.13.2 conda-forge::git=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgatk_0.0.19--py_0' :
        'quay.io/biocontainers/pypgatk:0.0.19--py_0' }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/git-lfs_1.5.2--0' :
        'quay.io/biocontainers/git-lfs:1.5.2--0' }"
    container "bitnami/git:2.30.0-debian-10-r27"

    when:
    params.cbioportal

    input:
    file cbioportal_config
    val cbioportal_study_id


    output:
    path('cbioportal_allstudies_data_mutations_mskcc.txt') ,emit: cbio_mutations
    path('cbioportal_allstudies_data_clinical_sample.txt') ,emit: cbio_samples

    script:
    if (cbioportal_study_id == "all")
        """
        git clone https://github.com/cBioPortal/datahub.git .
        git lfs install --local --skip-smudge
        git lfs pull -I public --include "data*clinical*sample.txt"
        git lfs pull -I public --include "data_mutations_mskcc.txt"
        cat public/*/data_mutations_mskcc.txt > cbioportal_allstudies_data_mutations_mskcc.txt
        cat public/*/*data*clinical*sample.txt | \\
            awk 'BEGIN{FS=OFS="\\t"}{if(\$1!~"#SAMPLE_ID"){gsub("#SAMPLE_ID", "\\nSAMPLE_ID");} print}' | \\
            awk 'BEGIN{FS=OFS="\\t"}{s=0; j=0; \\
                for(i=1;i<=NF;i++){ \\
                    if(\$i=="CANCER_TYPE_DETAILED") j=1; \\
                    if(\$i=="CANCER_TYPE") s=1; \\
                } \\
                if(j==1 && s==0){ \\
                    gsub("CANCER_TYPE_DETAILED", "CANCER_TYPE"); \\
                } \\
                print; \\
            }' \\
            > cbioportal_allstudies_data_clinical_sample.txt
        """
    else
        """
        pypgatk_cli.py cbioportal-downloader \\
            --config_file "$cbioportal_config" \\
            -d "$cbioportal_study_id"

        tar -xzvf database_cbioportal/${cbioportal_study_id}.tar.gz
        cat ${cbioportal_study_id}/data_mutations_mskcc.txt > cbioportal_allstudies_data_mutations_mskcc.txt
        cat ${cbioportal_study_id}/data_clinical_sample.txt | \\
            awk 'BEGIN{FS=OFS="\\t"}{if(\$1!~"#SAMPLE_ID"){gsub("#SAMPLE_ID", "\\nSAMPLE_ID");} print}' | \\
            awk 'BEGIN{FS=OFS="\\t"}{s=0; j=0; \\
            for(i=1;i<=NF;i++){ \\
                if(\$i=="CANCER_TYPE_DETAILED") j=1; if(\$i=="CANCER_TYPE") s=1; \\
            } \\
            if(j==1 && s==0){gsub("CANCER_TYPE_DETAILED", "CANCER_TYPE");} print;}' \\
            > cbioportal_allstudies_data_clinical_sample.txt
        """
}
