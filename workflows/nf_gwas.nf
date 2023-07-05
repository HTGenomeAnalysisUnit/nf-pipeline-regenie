//Set project_id
project_id = params.project

//Check required parameters for GWAS analysis
if (params.genotypes_imputed) {
  requiredParams = [
    'regenie_gwas_test',
    'regenie_min_imputation_score',
    'regenie_gwas_min_mac',
    'regenie_min_imputation_score'
  ]

  for (param in requiredParams) {
      if (params[param] == null || params[param] == '') {
        log.error "Parameter ${param} is required when running GWAS analysis."
        exit 1
      }
  }

  //Check specified regenie test is allowed
  def allowed_tests = ['additive', 'dominant', 'recessive']
  if (!(params.regenie_test in allowed_tests)){
    log.error "Test ${params.regenie_test} not supported. Allowed tests are $allowed_tests."
    exit 1
  }

  //Check imputed file format is allwed
  def allowed_input_formats = ['vcf', 'bgen', 'pgen', 'bed']
  if (!(params.genotypes_imputed_format in allowed_input_formats)){
    log.error "File format ${params.genotypes_imputed_format} not supported. Allowed formats are $allowed_input_formats."
    exit 1
  }
}

//Check required parameters for rare variants analysis
if (params.genotypes_rarevar) {
  requiredParams = [
    'rarevar_set_list_file',
    'rarevar_anno_file',
    'rarevar_mask_file',
    'regenie_rarevar_min_mac',
    'rarevars_aaf_bins',
    'rarevars_vc_test',
    'rarevars_vc_maxAAF',
    'regenie_build_mask'
  ]

  for (param in requiredParams) {
      if (params[param] == null || params[param] == '') {
        log.error "Parameter ${param} is required when running rare variant analysis."
        exit 1
      }
  }

  //Check imputed file format is allwed
  def allowed_input_formats = ['vcf', 'bgen', 'pgen', 'bed']
  if (!(params.genotypes_rarevar_format in allowed_input_formats)){
    log.error "File format ${params.genotypes_rarevar_format} not supported. Allowed formats are $allowed_input_formats."
    exit 1
  }
}

if (params.regenie_range != '' && ( params.step2_gwas_split || params.step2_rarevar_split )) {
  log.error "You cannot set regenie_range when step2_gwas_split and/or step2_rarevar_split is active"
  exit 1
}

//Set output and logs directories
if(params.outdir == null) {
  outdir = './'
} else {
  outdir = "${params.outdir}"
}

if (params.master_log_dir != null) {
  master_log_dir = "${params.master_log_dir}"
} else {
  master_log_dir = "${outdir}"
}

//Make chromosomes list
def chromosomes = []
params.chromosomes.split(',').each { element ->
    if (element.contains('-')) {
        def limits = element.split('-').collect { it.toInteger() }
        def range = limits[0] .. limits[1]
        chromosomes.addAll(range)
    } else {
        chromosomes << element
    }
}
chromosomes = chromosomes*.toString()

//Include modules
include { OPENING_LOG                 } from '../modules/local/opening_log' addParams(outdir: "$outdir")
include { REGENIE_STEP1_WF            } from '../subworkflow/regenie_step1' addParams(outdir: "$outdir", chromosomes: chromosomes)
include { REPORT_GWAS                 } from '../modules/local/report'  addParams(outdir: "$outdir")
include { REPORT_RAREVAR              } from '../modules/local/report'  addParams(outdir: "$outdir")

//gwas sub wf and modules
include { PREPARE_GENETIC_DATA as PREPARE_GWAS_DATA } from '../subworkflow/prepare_step2_data' addParams(outdir: "$outdir", genotypes_data: params.genotypes_imputed, input_format: params.genotypes_imputed_format, bgen_sample_file: params.imputed_sample_file, chromosomes: chromosomes, dosage_from: params.gwas_read_dosage_from)
include { SPLIT_GWAS_DATA_WF          } from '../subworkflow/split_data' addParams(outdir: "$outdir", chromosomes: chromosomes, input_format: params.genotypes_imputed_format)
include { PROCESS_GWAS_RESULTS_WF          } from '../subworkflow/process_results' addParams(outdir: "$outdir/results/gwas", logdir: "$outdir/logs", chromosomes: chromosomes, input_format: params.genotypes_imputed_format, rarevar_results: false, publish_filtered: false, annotation_min_log10p: params.annotation_min_log10p)
include { REGENIE_STEP2_WF as REGENIE_STEP2_GWAS_WF } from '../subworkflow/regenie_step2' addParams(outdir: "$outdir/results/gwas", logdir: "$outdir/logs", chromosomes: chromosomes, run_gwas: true, run_rarevar: false)

//rare variant sub wf and modules
include { PREPARE_GENETIC_DATA as PREPARE_RAREVARIANT_DATA } from '../subworkflow/prepare_step2_data' addParams(outdir: "$outdir", genotypes_data: params.genotypes_rarevar, input_format: params.genotypes_rarevar_format, bgen_sample_file: params.rarevar_sample_file, chromosomes: chromosomes, dosage_from: params.rarevar_read_dosage_from)
include { SPLIT_RAREVARIANT_DATA_WF        } from '../subworkflow/split_data' addParams(outdir: "$outdir", chromosomes: chromosomes, input_format: params.genotypes_imputed_format)
include { PROCESS_RAREVAR_RESULTS_WF        } from '../subworkflow/process_results' addParams(outdir: "$outdir/results/rarevar", logdir: "$outdir/logs", chromosomes: chromosomes, rarevar_results: true, publish_filtered: true, annotation_min_log10p: params.rarevar_min_log10p)
include { REGENIE_STEP2_WF as REGENIE_STEP2_RAREVAR_WF } from '../subworkflow/regenie_step2' addParams(outdir: "$outdir/results/rarevar", logdir: "$outdir/logs", chromosomes: chromosomes, run_gwas: false, run_rarevar: true)

//==== WORKFLOW ====
workflow NF_GWAS {   
  take:
  project_data //[project_id, pheno_file, pheno_meta(cols, binary, model), covar_file, covar_meta(cols, cat_cols)]
  input_validation_logs //[project_id, pheno_validation_log, covar_validation_log]

  main:
  //==== SET WORKFLOW runName ====
  workflow.runName = "${workflow.runName}-${project_data[0]}"

  //==== OPENING LOG ====
  OPENING_LOG(project_data)

  //TODO: Combine genotyped_plink_ch with project data and perform step 1 by project
  //TODO: results of step 1 are mapped into project_data
  //==== STEP 1 ====
  //Set input channel for step 1
  genotyped_files = Channel.fromFilePairs("${params.genotypes_array}.{bed,bim,fam}", size: 3, flat: true)
    .combine(project_data.map { it[0] } )
    .set {genotyped_plink_ch}

  REGENIE_STEP1_WF (
    genotyped_plink_ch,
    project_data
  )
  

  //==== STEP2 AND REPORTS - GWAS ====
  if (params.genotypes_imputed) {
    //Prepare data for step 2
    PREPARE_GWAS_DATA()
    if (params.step2_gwas_split) {
      SPLIT_GWAS_DATA_WF(PREPARE_GWAS_DATA.out.processed_genotypes)
      gwas_data_input_ch = SPLIT_GWAS_DATA_WF.out.processed_genotypes
    } else {
      gwas_data_input_ch = PREPARE_GWAS_DATA.out.processed_genotypes
        .map { tuple (it[0], it[1], it[2], it[3], it[4], "SINGLE_CHUNK") }
    }

    //TODO: Combine gwas_data_input_ch with step1 project data. Subsequent steps runs by project data.
    //Run regenie step 2
    REGENIE_STEP2_GWAS_WF(
      gwas_data_input_ch,
      REGENIE_STEP1_WF.out.regenie_step1_out,
      VALIDATE_PHENOTYPES.out.phenotypes_file_validated,
      covariates_file_validated,
      CACHE_JBANG_SCRIPTS.out.regenie_log_parser_jar
    )

    //Process results - filtering and annotation
    PROCESS_GWAS_RESULTS_WF(REGENIE_STEP2_GWAS_WF.out.regenie_results, PREPARE_GWAS_DATA.out.processed_genotypes)

    //==== GENERATE HTML REPORTS ====
    if (params.make_report) {
      gwas_report_template = file("$projectDir/reports/gwas_report_template.Rmd",checkIfExists: true)
      REPORT_GWAS (
          PROCESS_GWAS_RESULTS_WF.out.processed_results,
          VALIDATE_PHENOTYPES.out.phenotypes_file_validated,
          gwas_report_template,
          VALIDATE_PHENOTYPES.out.phenotypes_file_validated_log,
          covariates_file_validated_log,
          REGENIE_STEP1_WF.out.regenie_step1_parsed_logs.collect(),
          REGENIE_STEP2_GWAS_WF.out.regenie_log
      )
    }
  }
  
  //==== STEP2 AND REPORTS - RARE VARIANTS ====
  if (params.genotypes_rarevar) {
    //Prpare data for step 2
    PREPARE_RAREVARIANT_DATA()
    
    if (params.step2_rarevar_split) {
      SPLIT_RAREVARIANT_DATA_WF(PREPARE_RAREVARIANT_DATA.out.processed_genotypes)
      rare_variants_input_ch = SPLIT_RAREVARIANT_DATA_WF.out.processed_genotypes
    } else {
      rare_variants_input_ch = PREPARE_RAREVARIANT_DATA.out.processed_genotypes
        .map { tuple (it[0], it[1], it[2], it[3], it[4], "SINGLE_CHUNK") }
    }

    //TODO: Combine rare_variants_input_ch with step1 project data. Subsequent steps runs by project data.
    //Run regenie step 2
    REGENIE_STEP2_RAREVAR_WF(
      rare_variants_input_ch,
      REGENIE_STEP1_WF.out.regenie_step1_out,
      VALIDATE_PHENOTYPES.out.phenotypes_file_validated,
      covariates_file_validated,
      CACHE_JBANG_SCRIPTS.out.regenie_log_parser_jar
    )

    //Process results - filtering
    PROCESS_RAREVAR_RESULTS_WF(REGENIE_STEP2_RAREVAR_WF.out.regenie_results)

    //==== GENERATE HTML REPORTS ====
    if (params.make_report) {
      rarevar_report_template = file("$projectDir/reports/rare_vars_report_template.Rmd",checkIfExists: true)
      REPORT_RAREVAR (
        PROCESS_RAREVAR_RESULTS_WF.out.processed_results,
        VALIDATE_PHENOTYPES.out.phenotypes_file_validated,
        rarevar_report_template,
        VALIDATE_PHENOTYPES.out.phenotypes_file_validated_log,
        covariates_file_validated_log,
        REGENIE_STEP1_WF.out.regenie_step1_parsed_logs.collect(),
        REGENIE_STEP2_RAREVAR_WF.out.regenie_log
      )
    }
  }
}

workflow.onComplete {
  //==== SAVE CONFIGURATION ====
  pipeline_log_dir = file("$outdir/analysis_config")
  pipeline_log_dir.mkdirs()
  
  def msg="""\
  ## Pipeline information ##
  Pipeline version = ${workflow.manifest.version}
  Pipeline repo = ${workflow.manifest.homePage}
  Nextflow version required = ${workflow.manifest.nextflowVersion}

  ## Execution summary ##
  Execution status: ${ workflow.success ? 'OK' : 'failed' }
  Execution name: ${workflow.runName}
  Execution ID: ${workflow.sessionId}
  Execution date: ${workflow.start}
  Executed by: ${workflow.userName}
  Commandline: ${workflow.commandLine} 
  Container: ${workflow.container}
  """
  .stripIndent()

  params_log = file("${pipeline_log_dir}/pipeline_execution.log")
  params_log.append(msg)

  for (f in workflow.configFiles) {
    config_file = file(f)
    f.copyTo("$pipeline_log_dir")
  }

  nextflow_log = file("${workflow.launchDir}/.nextflow.log")
  nextflow_log.copyTo("${pipeline_log_dir}/nextflow.log")

  versions_yml = file("${projectDir}/conf/main_tools_versions.yml")
  versions_yml.copyTo("${pipeline_log_dir}/main_tools_versions.yml")

  //if (params.master_log_dir != null) {
    master_log = file("${master_log_dir}/job_execution_summary.log")
    def master_log_msg="${params.project}\t${workflow.launchDir}\t${ workflow.success ? 'OK' : 'failed' }\n"
    master_log.append(master_log_msg)
  //}

  //CLOSE MESSAGE
  println "Results location: ${ outdir }"
}

workflow.onError {
  if (params.master_log_dir != null) {
    error_log = file("${master_log_dir}/${project_id}_error.log")
    def error_log_msg="""\
    ### ERROR MESSAGE ###
    ${workflow.errorMessage}

    ### ERROR REPORT ###
    ${workflow.errorReport}
    """
    .stripIndent()
    error_log.append(error_log_msg)

    master_log = file("${master_log_dir}/job_execution_summary.log")
    def master_log_msg="${params.project}\t${workflow.launchDir}\t${ workflow.success ? 'OK' : 'failed' }\n"
    master_log.append(master_log_msg)
  }
}
