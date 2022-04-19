#!/usr/bin/env nextflow
/*
========================================================================================
    genepi/nf-gwas
========================================================================================
    GitLab : https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie
    Author : Edoardo Giacopuzzi
    Based on github : https://github.com/genepi/nf-gwas
    ---------------------------
*/

nextflow.enable.dsl = 2

include { NF_GWAS } from './workflows/nf_gwas'

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

workflow {
    NF_GWAS ()
}
