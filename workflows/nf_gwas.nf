//Check general required parameters
requiredParams = [
  'project', 'genotypes_array', 'genotypes_build',
  'phenotypes_filename', 'phenotypes_columns', 'phenotypes_binary_trait',
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
        exit 1, "Parameter ${param} is required when running GWAS analysis."
      }
  }

  //Check specified regenie test is allowed
  def allowed_tests = ['additive', 'dominant', 'recessive']
  if (!(params.regenie_test in allowed_tests)){
    exit 1, "Test ${params.regenie_test} not supported. Allowed tests are $allowed_tests."
  }

  //Check imputed file format is allwed
  def allowed_input_formats = ['vcf', 'bgen', 'pgen', 'bed']
  if (!(params.genotypes_imputed_format in allowed_input_formats)){
    exit 1, "File format ${params.genotypes_imputed_format} not supported. Allowed formats are $allowed_input_formats."
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
        exit 1, "Parameter ${param} is required when running rare variant analysis."
      }
  }

  //Check imputed file format is allwed
  def allowed_input_formats = ['vcf', 'bgen', 'pgen', 'bed']
  if (!(params.genotypes_rarevar_format in allowed_input_formats)){
    exit 1, "File format ${params.genotypes_rarevar_format} not supported. Allowed formats are $allowed_input_formats."
  }
}

if (params.regenie_range != '' && ( params.step2_gwas_split || params.step2_rarevar_split )) {
  exit 1, "You cannot set regenie_range when step2_gwas_split and/or step2_rarevar_split is active"
}

//Set output and logs directories
if(params.outdir == null) {
  outdir = "${params.project}"
} else {
  outdir = "${params.outdir}/${params.project}"
}

if (params.master_log_dir != null) {
  master_log_dir = "${params.master_log_dir}"
} else {
  master_log_dir = "${outdir}"
}

project_id = params.project
phenotypes_array = params.phenotypes_columns.trim().split(',')

covariates_array= []
if(!params.covariates_columns.isEmpty()){
  covariates_array = params.covariates_columns.trim().split(',')
}

//Check accessory scripts
regenie_log_parser_java  = file("$projectDir/bin/RegenieLogParser.java", checkIfExists: true)
regenie_filter_java = file("$projectDir/bin/RegenieFilter.java", checkIfExists: true)
regenie_validate_input_java = file("$projectDir/bin/RegenieValidateInput.java", checkIfExists: true)

//Phenotypes
phenotypes_file = file(params.phenotypes_filename, checkIfExists: true)
phenotypes = Channel.from(phenotypes_array)

//Covariates
if (params.covariates_filename == 'NO_COV_FILE') {
  covar_tmp_file = file("${workflow.workDir}/NO_COV_FILE")
  covar_tmp_file.append('')
  covariates_file = file("$covar_tmp_file")
} else {
  covariates_file = file(params.covariates_filename)
}
if (params.covariates_filename != 'NO_COV_FILE' && !covariates_file.exists()){
  exit 1, "Covariate file ${params.covariates_filename} not found."
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

//Set input channel for step 1
Channel.fromFilePairs("${params.genotypes_array}.{bim,bed,fam}", size: 3).set {genotyped_plink_ch}

//Include modules
include { CACHE_JBANG_SCRIPTS         } from '../modules/local/cache_jbang_scripts'
include { VALIDATE_PHENOTYPES         } from '../modules/local/validate_phenotypes' addParams(outdir: "$outdir")
include { VALIDATE_COVARIATS          } from '../modules/local/validate_covariates' addParams(outdir: "$outdir")
include { REGENIE_STEP1_WF            } from '../subworkflow/regenie_step1' addParams(outdir: "$outdir", chromosomes: chromosomes)
include { REPORT_GWAS                 } from '../modules/local/report'  addParams(outdir: "$outdir")
include { REPORT_RAREVAR              } from '../modules/local/report'  addParams(outdir: "$outdir")


//if (params.run_gwas) {
  include { PREPARE_GENETIC_DATA as PREPARE_GWAS_DATA } from '../subworkflow/prepare_step2_data' addParams(outdir: "$outdir", genotypes_data: params.genotypes_imputed, input_format: params.genotypes_imputed_format, bgen_sample_file: params.imputed_sample_file, chromosomes: chromosomes, dosage_from: params.gwas_read_dosage_from)
  include { SPLIT_GWAS_DATA_WF          } from '../subworkflow/split_data' addParams(outdir: "$outdir", chromosomes: chromosomes, input_format: params.genotypes_imputed_format)
  include { PROCESS_GWAS_RESULTS_WF          } from '../subworkflow/process_results' addParams(outdir: "$outdir/results/gwas", logdir: "$outdir/logs", chromosomes: chromosomes, input_format: params.genotypes_imputed_format, rarevar_results: false, publish_filtered: false)
  include { REGENIE_STEP2_WF as REGENIE_STEP2_GWAS_WF } from '../subworkflow/regenie_step2' addParams(outdir: "$outdir/results/gwas", logdir: "$outdir/logs", chromosomes: chromosomes, run_gwas: true, run_rarevar: false)
//}
//if (params.run_rare_variants) {
  include { PREPARE_GENETIC_DATA as PREPARE_RAREVARIANT_DATA } from '../subworkflow/prepare_step2_data' addParams(outdir: "$outdir", genotypes_data: params.genotypes_rarevar, input_format: params.genotypes_rarevar_format, bgen_sample_file: params.rarevar_sample_file, chromosomes: chromosomes, dosage_from: params.rarevar_read_dosage_from)
  include { SPLIT_RAREVARIANT_DATA_WF        } from '../subworkflow/split_data' addParams(outdir: "$outdir", chromosomes: chromosomes, input_format: params.genotypes_imputed_format)
  include { PROCESS_RAREVAR_RESULTS_WF        } from '../subworkflow/process_results' addParams(outdir: "$outdir/results/rarevar", logdir: "$outdir/logs", chromosomes: chromosomes, rarevar_results: true, publish_filtered: true)
  include { REGENIE_STEP2_WF as REGENIE_STEP2_RAREVAR_WF } from '../subworkflow/regenie_step2' addParams(outdir: "$outdir/results/rarevar", logdir: "$outdir/logs", chromosomes: chromosomes, run_gwas: false, run_rarevar: true)
//}


//==== WORKFLOW ====
workflow NF_GWAS {   

  //==== INITIAL LOGGING OF PARAMETERS ====
  log_params = [ 
  'genotypes_array',
  'genotypes_imputed', 'genotypes_imputed_format',
  'genotypes_rarevar', 'genotypes_rarevar_format',
  'genotypes_build',
  'phenotypes_filename','phenotypes_columns','phenotypes_binary_trait',
  'covariates_filename','covariates_columns',
  'chromosomes',
  'regenie_test','annotation_min_log10p',
  'save_step1_predictions'
  ]
  log_params_string = []
  for (p in log_params) {
      log_params_string.add("$p : " + params[p])
  }

log.info"""\
==========================================================
  REGENIE GWAS - PROJECT ${params.project} - NF PIPELINE    
==========================================================
PROJECT ID      : ${params.project}
OUTDIR          : ${outdir}

${log_params_string.join('\n')}
==========================================================
Please report issues to:
https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie
or contact: edoardo.giacopuzzi@fht.org
"""

  //==== VALIDATE COVARS AND PHENO TABLES ====
  CACHE_JBANG_SCRIPTS (
      regenie_log_parser_java,
      regenie_filter_java,
      regenie_validate_input_java
  )

  VALIDATE_PHENOTYPES (
      phenotypes_file,
      CACHE_JBANG_SCRIPTS.out.regenie_validate_input_jar
  )

  if(covariates_file.exists() && covariates_file.name != 'NO_COV_FILE') {
      VALIDATE_COVARIATS (
        covariates_file,
        CACHE_JBANG_SCRIPTS.out.regenie_validate_input_jar
      )

      covariates_file_validated = VALIDATE_COVARIATS.out.covariates_file_validated
      covariates_file_validated_log = VALIDATE_COVARIATS.out.covariates_file_validated_log

  } else {

    // set covariates_file to default value
    covariates_file_validated = covariates_file
    covariates_file_validated_log = Channel.fromPath("NO_COV_LOG")
  }

  //==== STEP 1 ====
  REGENIE_STEP1_WF (
    genotyped_plink_ch,
    VALIDATE_PHENOTYPES.out.phenotypes_file_validated,
    covariates_file_validated,
    CACHE_JBANG_SCRIPTS.out.regenie_log_parser_jar
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
