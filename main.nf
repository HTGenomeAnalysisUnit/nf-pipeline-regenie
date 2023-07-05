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

//Check general required parameters
requiredParams = [
  'project', 'genotypes_array', 'genotypes_build',
  'chromosomes',
  'prune_enabled',
  'prune_maf',
  'prune_window_kbsize',
  'prune_step_size',
  'prune_r2_threshold',
  'qc_maf',
  'qc_mac',
  'qc_geno',
  'qc_hwe',
  'qc_mind',
  'regenie_bsize_step1',
  'step1_n_chunks',
  'regenie_bsize_step2',
  'annotation_min_log10p',
  'annotation_interval_kb'
]

for (param in requiredParams) {
    if (params[param] == null || params[param] == '') {
      exit 1, "Parameter ${param} is required."
    }
}

include { PREPARE_PROJECT } from './workflows/prepare_project'
include { NF_GWAS } from './workflows/nf_gwas'

workflow {
    //==== INITIAL LOGGING OF PARAMETERS ====
    log_params = [ 
    'genotypes_array',
    'genotypes_imputed', 'genotypes_imputed_format',
    'genotypes_rarevar', 'genotypes_rarevar_format',
    'genotypes_build',
    'chromosomes',
    'annotation_min_log10p',
    'save_step1_predictions'
    ]
    log_params_string = []
    for (p in log_params) {
    log_params_string.add("$p : " + params[p])
    }

log.info"""\
==========================================================
  REGENIE GWAS - SHARED PARAMETERS - NF PIPELINE    
==========================================================

${global_parameters.join('\n')}
==========================================================
Please report issues to:
https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie
or contact: edoardo.giacopuzzi@fht.org
"""
    //==== PREPARE PROJECT INPUTS ====
    PREPARE_PROJECT()

    //==== RUN NF-GWAS ====
    // project_data = [project_id, pheno_file, pheno_meta(cols, binary, model), covar_file, covar_meta(cols, cat_cols)]
    NF_GWAS(PREAPARE_PROJECT.out.project_data, PREPARE_PROJECT.out.input_validation_logs)
    
}
