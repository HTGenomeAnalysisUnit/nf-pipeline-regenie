#!/usr/bin/env nextflow
/*
========================================================================================
    nf-pipeline-regenie
========================================================================================
    GitHub : https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie
    Author : Edoardo Giacopuzzi
    ---------------------------
*/

nextflow.enable.dsl = 2
import java.text.SimpleDateFormat

/*
======================================================================
    INITIALIZATION
======================================================================
*/

//Check general required parameters
requiredParams = [
  'project', 'genotypes_build',
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

if (!(params.regenie_skip_predictions || params.regenie_premade_predictions)) {
  if (params.genotypes_array == null || params.genotypes_array == '') {
    exit 1, "Parameter genotypes_array is required when regenie_skip_predictions or regenie_premade_predictions are not set"
  }
}

if (params.regenie_range != '' && ( params.step2_gwas_split || params.step2_rarevar_split )) {
  log.error "You cannot set regenie_range when step2_gwas_split and/or step2_rarevar_split is active"
  exit 1
}

//Set output and logs directories
def date = new Date()
def formatted_date = new SimpleDateFormat("yyyyMMdd_HHmmss")
if(params.outdir == null) {
  outdir = "pipeline_results_${formatted_date}"
} else {
  outdir = "${params.outdir}"
}

if (params.master_log_dir == null) {
  master_log_dir = "${outdir}"
} else {
  master_log_dir = "${params.master_log_dir}"
}

include { PREPARE_PROJECT       } from './workflows/prepare_project'  addParams(outdir: outdir, logdir: master_log_dir)
include { RUN_VARIANT_ANALYSIS  } from './workflows/variant_analysis' addParams(outdir: outdir, logdir: master_log_dir)

/*
======================================================================
    WORKFLOW
======================================================================
*/

workflow {
  //==== SET WORKFLOW runName ====
  workflow.runName = "${params.project}-${workflow.runName}"

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
  
  global_parameters = []
  for (p in log_params) {
    global_parameters.add("$p : " + params[p])
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
  RUN_VARIANT_ANALYSIS(PREPARE_PROJECT.out.project_data, PREPARE_PROJECT.out.input_validation_logs)
    
}
