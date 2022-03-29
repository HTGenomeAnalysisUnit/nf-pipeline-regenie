#!/usr/bin/env nextflow
/*
========================================================================================
    genepi/nf-gwas
========================================================================================
    Github : https://github.com/genepi/nf-gwas
    Author: Sebastian Schönherr & Lukas Forer
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
    log.info "PARALLEL CHROM VERSION"
    NF_GWAS ()
}
