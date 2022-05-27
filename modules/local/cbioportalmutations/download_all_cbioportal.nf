/**
 * Download all cBioPortal studies using git-lfs
*/
process DOWNLOAD_ALL_CBIOPORTAL {

    when:
    params.cbioportal

    output:
    path('cbioportal_allstudies_data_mutations_mskcc.txt') ,emit: cbio_mutations
    path('cbioportal_allstudies_data_clinical_sample.txt') ,emit: cbio_samples

    script:
    if (params.cbioportal_study_id == "all")
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
            -d "$params.cbioportal_study_id"

        tar -xzvf database_cbioportal/${params.cbioportal_study_id}.tar.gz
        cat ${params.cbioportal_study_id}/data_mutations_mskcc.txt > cbioportal_allstudies_data_mutations_mskcc.txt
        cat ${params.cbioportal_study_id}/data_clinical_sample.txt | \\
            awk 'BEGIN{FS=OFS="\\t"}{if(\$1!~"#SAMPLE_ID"){gsub("#SAMPLE_ID", "\\nSAMPLE_ID");} print}' | \\
            awk 'BEGIN{FS=OFS="\\t"}{s=0; j=0; \\
            for(i=1;i<=NF;i++){ \\
                if(\$i=="CANCER_TYPE_DETAILED") j=1; if(\$i=="CANCER_TYPE") s=1; \\
            } \\
            if(j==1 && s==0){gsub("CANCER_TYPE_DETAILED", "CANCER_TYPE");} print;}' \\
            > cbioportal_allstudies_data_clinical_sample.txt
        """
}