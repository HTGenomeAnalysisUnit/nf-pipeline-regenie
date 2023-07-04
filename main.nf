#!/usr/bin/env nextflow
/*
========================================================================================
    nf-pipeline-regenie
========================================================================================
    GitHub : https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie
    Author : Edoardo Giacopuzzi
    Original concept GitHub : https://github.com/genepi/nf-gwas
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
