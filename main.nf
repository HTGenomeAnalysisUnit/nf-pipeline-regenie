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

if (params.with_master) {
    include { SETUP_MULTIPLE_RUNS } from './workflows/setup_runs'
} else {
    include { NF_GWAS } from './workflows/nf_gwas'
}

workflow {
    if (params.with_master) {
        println "Using a master mode"
        SETUP_MULTIPLE_RUNS()
    } else {
        NF_GWAS()
    }
}
