//Set variables
def allowed_input_formats = ['vcf', 'bgen', 'pgen', 'bed']
def allowed_rarevar_stats = ["FDR_bygroup","FDR_alltests","BONF_bygroup","BONF_alltests"]
quarto_report_css = file("$projectDir/reports/quarto_report.css", checkIfExists: true)
make_chunks_sql = file("$projectDir/bin/get_intervals.sql", checkIfExists: true)

//Check required parameters for GWAS analysis
if (params.genotypes_imputed) {
  requiredParams = [
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

  //Check imputed file format is allowed
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

  //Check imputed file format is allowed
  if (!(params.genotypes_rarevar_format in allowed_input_formats)){
    log.error "File format ${params.genotypes_rarevar_format} not supported. Allowed formats are $allowed_input_formats."
    exit 1
  }

  //Check rare variant stats for report is allowed
  if (!(params.rarevar_stat_test in allowed_rarevar_stats)){
    log.error "Rare variants stat column ${params.rarevar_stat_test} not supported. Allowed formats are $allowed_rarevar_stats."
    exit 1
  }
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
include { REGENIE_STEP1_WF            } from '../subworkflow/regenie_step1' addParams(chromosomes: chromosomes)
include { REPORT_GWAS                 } from '../modules/local/report'
include { REPORT_RAREVAR              } from '../modules/local/report'

//gwas sub wf and modules
include { PREPARE_GENETIC_DATA as PREPARE_GWAS_DATA } from '../subworkflow/prepare_step2_data' addParams(genotypes_data: params.genotypes_imputed, input_format: params.genotypes_imputed_format, bgen_sample_file: params.imputed_sample_file, chromosomes: chromosomes, dosage_from: params.gwas_read_dosage_from)
include { SPLIT_GWAS_DATA_WF          } from '../subworkflow/split_data' addParams(chromosomes: chromosomes, input_format: params.genotypes_imputed_format)
include { PROCESS_GWAS_RESULTS_WF          } from '../subworkflow/process_results' addParams(chromosomes: chromosomes, rarevar_results: false, input_format: params.genotypes_imputed_format, input_files: params.genotypes_imputed, annotation_min_log10p: params.annotation_min_log10p)
include { REGENIE_STEP2_WF as REGENIE_STEP2_GWAS_WF } from '../subworkflow/regenie_step2' addParams(chromosomes: chromosomes, run_gwas: true, run_rarevar: false)

//rare variant sub wf and modules
include { PREPARE_GENETIC_DATA as PREPARE_RAREVARIANT_DATA } from '../subworkflow/prepare_step2_data' addParams(genotypes_data: params.genotypes_rarevar, input_format: params.genotypes_rarevar_format, bgen_sample_file: params.rarevar_sample_file, chromosomes: chromosomes, dosage_from: params.rarevar_read_dosage_from)
include { SPLIT_RAREVARIANT_DATA_WF        } from '../subworkflow/split_data' addParams(chromosomes: chromosomes, input_format: params.genotypes_imputed_format)
include { PROCESS_RAREVAR_RESULTS_WF        } from '../subworkflow/process_results' addParams(chromosomes: chromosomes, rarevar_results: true, annotation_min_log10p: params.rarevar_min_log10p)
include { REGENIE_STEP2_WF as REGENIE_STEP2_RAREVAR_WF } from '../subworkflow/regenie_step2' addParams(chromosomes: chromosomes, run_gwas: false, run_rarevar: true)

//==== WORKFLOW ====
workflow RUN_VARIANT_ANALYSIS {   
  take:
  project_data //[project_id, pheno_file, pheno_meta(cols, binary, model), covar_file, covar_meta(cols, cat_cols), [accessory_files]]
  input_validation_logs //[project_id, pheno_validation_log, covar_validation_log]

  main:
  //==== STEP 1 ====
  //Set input channel for step 1
  genotyped_files = Channel.fromFilePairs("${params.genotypes_array}.{bed,bim,fam}", size: 3, flat: true)
    .map{ tuple(it[1], it[2], it[3])}
  genotyped_plink_ch = project_data.map { it[0] }.combine(genotyped_files)
  
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
      processed_gwas_data_ch = SPLIT_GWAS_DATA_WF.out.processed_genotypes
    } else {
      //If there is no split, we need to count the number of files provided in input
      //In this way, if the user provide a dataset spread across multiple files we can merge results correctly with groupTuple
      count_chunks = PREPARE_GWAS_DATA.out.count()
      processed_gwas_data_ch = PREPARE_GWAS_DATA.out.processed_genotypes
        .map { tuple (it[0], it[1], it[2], it[3], it[4], "SINGLE_CHUNK") }
        .combine(count_chunks)
    }

    //Run regenie step 2
    gwas_data_input_ch = project_data
      .join(REGENIE_STEP1_WF.out.regenie_step1_out)
      .combine(processed_gwas_data_ch)

    REGENIE_STEP2_GWAS_WF( gwas_data_input_ch )

    //Process results - filtering and annotation
    PROCESS_GWAS_RESULTS_WF(REGENIE_STEP2_GWAS_WF.out.regenie_results, PREPARE_GWAS_DATA.out.processed_genotypes)

    //==== GENERATE HTML REPORTS ====
    logs_ch = project_data.map{ tuple(it[0], it[1]) }
      .join(input_validation_logs)
      .join(REGENIE_STEP1_WF.out.regenie_step1_parsed_logs)
      .join(REGENIE_STEP2_GWAS_WF.out.regenie_log)

    report_input_ch = logs_ch.combine(PROCESS_GWAS_RESULTS_WF.out.processed_results, by:0)
      //[val(project_id), path(phenotype_file), path(phenotype_log), path(covariate_log), path(step1_log), path(step2_log), val(phenotype), path(regenie_merged_results), path(annotated_tophits), path(annotated_toploci)]

    if (params.make_report) {
      gwas_report_template = file("$projectDir/reports/gwas_report_template.qmd", checkIfExists: true)
      REPORT_GWAS (
          report_input_ch,
          gwas_report_template,
          quarto_report_css
      )
    }
  }
  
  //==== STEP2 AND REPORTS - RARE VARIANTS ====
  if (params.genotypes_rarevar) {
    //Prpare data for step 2
    PREPARE_RAREVARIANT_DATA()
    
    if (params.step2_rarevar_split) {
      SPLIT_RAREVARIANT_DATA_WF(PREPARE_RAREVARIANT_DATA.out.processed_genotypes)
      processed_rarevar_data_ch = SPLIT_RAREVARIANT_DATA_WF.out.processed_genotypes
    } else {
      processed_rarevar_data_ch = PREPARE_RAREVARIANT_DATA.out.processed_genotypes
        .map { tuple (it[0], it[1], it[2], it[3], it[4], "SINGLE_CHUNK") }
    }

    rarevar_data_input_ch = project_data
      .join(REGENIE_STEP1_WF.out.regenie_step1_out)
      .combine(processed_rarevar_data_ch)

    //Run regenie step 2
    REGENIE_STEP2_RAREVAR_WF( rarevar_data_input_ch )

    //Process results - filtering
    PROCESS_RAREVAR_RESULTS_WF(REGENIE_STEP2_RAREVAR_WF.out.regenie_results)

    //==== GENERATE HTML REPORTS ====
    report_input_ch = project_data.map{ tuple(it[0], it[1], it[4]) }
      .join(input_validation_logs)
      .join(REGENIE_STEP1_WF.out.regenie_step1_parsed_logs)
      .join(REGENIE_STEP2_RAREVAR_WF.out.regenie_log)
      .combine(PROCESS_RAREVAR_RESULTS_WF.out.processed_results, by:0)

    if (params.make_report) {
      rarevar_report_template = file("$projectDir/reports/rare_vars_report_template.qmd",checkIfExists: true)
      REPORT_RAREVAR (
          report_input_ch,
          rarevar_report_template,
          quarto_report_css
      )
    }
  }
}

workflow.onComplete {
  //==== SAVE CONFIGURATION ====
  pipeline_log_dir = file("${params.logdir}/analysis_config")
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

  //CLOSE MESSAGE
  log.info """
  ==============================
  ANALYSIS COMPLETE!
  Execution status: ${ workflow.success ? 'OK' : 'FAILED' }
  Results location: ${ params.outdir }
  """.stripIndent()
}