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

nextflow.enable.dsl = 2

//WorkflowMain.initialise(workflow, params, log)

 include { PGDB } from './workflows/pgdb'

//
// WORKFLOW: Run main nf-core/pgdb analysis pipeline
//
 workflow {
    PGDB()   
}

/*
========================================================================================
    THE END
========================================================================================
*/