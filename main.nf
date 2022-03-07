#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/snrnaseq
========================================================================================
    Github : https://github.com/nf-core/snrnaseq
    Website: https://nf-co.re/snrnaseq
    Slack  : https://nfcore.slack.com/channels/snrnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { SNRNASEQ } from './workflows/snrnaseq'

//
// WORKFLOW: Run main nf-core/snrnaseq analysis pipeline
//
workflow NFCORE_SNRNASEQ {
    SNRNASEQ ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_SNRNASEQ ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
