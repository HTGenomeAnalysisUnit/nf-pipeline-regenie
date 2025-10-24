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
    INITIALIZATION AND VALIDATION
======================================================================
*/

// Include nf-schema plugin functions
include { paramsSummaryLog; paramsSummaryMap; paramsHelp } from 'plugin/nf-schema'

/*
======================================================================
    PARAMETER VALIDATION AND HELP
======================================================================
*/

// Show help message if requested
if (params.help) {
    def String command = "nextflow run ${manifest.name} --project <PROJECT_ID> --genotypes_build <BUILD>"
    log.info paramsHelp(command, parameters_schema: "$projectDir/nextflow_schema.json")
    exit 0
}

// Print parameter summary log with enhanced formatting
def summary_params = paramsSummaryMap(workflow, parameters_schema: "$projectDir/nextflow_schema.json")
log.info paramsSummaryLog(workflow, parameters_schema: "$projectDir/nextflow_schema.json")

// Log key pipeline information
log.info """\
==========================================================
PIPELINE INFORMATION
==========================================================
Pipeline : ${manifest.name}
Version  : ${manifest.version} 
Git info : ${workflow.repository} - ${workflow.revision} [${workflow.commitId}]
Command  : ${workflow.commandLine}
Profile  : ${workflow.profile}
Work dir : ${workflow.workDir}
==========================================================
"""

// Validate parameters using nf-schema - this will handle required parameter checks
// Note: The manual validation has been replaced by the schema validation which is more robust

//Additional conditional validation that cannot be easily expressed in JSON schema
if (!(params.regenie_skip_predictions || params.regenie_premade_predictions)) {
  if (params.genotypes_array == null || params.genotypes_array == '') {
    exit 1, "Parameter genotypes_array is required when regenie_skip_predictions or regenie_premade_predictions are not set"
  }
}

if ((params.regenie_range != '' || params.regenie_extract_snps != '') && params.step2_gwas_split ) {
  log.error "You cannot set regenie_range or regenie_extract_snps when step2_gwas_split is active"
  exit 1
}

if ((params.regenie_range != '' || params.regenie_extract_genes != '') && params.step2_rarevar_split ) {
  log.error "You cannot set regenie_range or regenie_extract_genes when step2_rarevar_split is active"
  exit 1
}

//Set output and logs directories
if(params.outdir == null) {
  outdir = "${params.project}_output"
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

  //==== ENHANCED PARAMETER LOGGING ====
  // The comprehensive parameter summary is now provided by nf-schema at the beginning
  log.info """\
==========================================================
  REGENIE GWAS - HIGH-SPEED PIPELINE    
==========================================================
Pipeline: ${manifest.name} v${manifest.version}
==========================================================
Please report issues to:
https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie
or contact: edoardo.giacopuzzi@fht.org
==========================================================
"""

  //==== PREPARE PROJECT INPUTS ====
  PREPARE_PROJECT()

  //==== RUN NF-GWAS ====
  // project_data = [project_id, pheno_file, pheno_meta(cols, binary, model), covar_file, covar_meta(cols, cat_cols)]
  RUN_VARIANT_ANALYSIS(PREPARE_PROJECT.out.project_data, PREPARE_PROJECT.out.input_validation_logs)
    
}
