pheno_chunker_r = file("$projectDir/bin/pheno_chunker.R", checkIfExists: true)
prepare_projects_py = file("$projectDir/bin/prepare_projects.py", checkIfExists: true)
conf_template = file("$projectDir/templates/run_parameters.conf", checkIfExists: true)

if(!params.master_outdir) {
    exit 1, "Please provide master output directory setting master_outdir parameter"
}
master_outdir = file(params.master_outdir)
master_outdir.mkdirs()

shared_config_file = file(params.shared_config_file, checkIfExists: true)
models_table = file(params.models_table, checkIfExists: true)
traits_table = file(params.traits_table, checkIfExists: true)

models_table.copyTo("$master_outdir/models.tsv")

if (params.master_log_dir == null) {
    master_log_dir = file("${master_outdir}/execution_log")
} else {
    master_log_dir = file(params.master_log_dir)
}
master_log_dir.mkdirs()

for (x in ["fam", "bim", "bed"]) {
    file("${params.genotypes_array}.${x}", checkIfExists: true)
}

fam_file = file("${params.genotypes_array}.fam", checkIfExists: true)

include { PREPARE_PROJECTS } from '../modules/local/prepare_projects' addParams(outdir: "$master_outdir")
include { SUBMIT_GWAS_RUN } from '../modules/local/submit_gwas_run' addParams(outdir: "$master_outdir", master_log_dir: "$master_log_dir")

workflow SETUP_MULTIPLE_RUNS {
    PREPARE_PROJECTS(pheno_chunker_r, prepare_projects_py, traits_table, models_table, fam_file, conf_template) // This runs phenotyper and prepares the project configuration
    SUBMIT_GWAS_RUN(PREPARE_PROJECTS.out.chunks.flatten(), shared_config_file) //This spawn a nextflow process for each configured run. Probably better 1 or 2 at a time
}

workflow.onComplete {
  //==== SAVE CONFIGURATION ====
  pipeline_log_dir = file("$master_outdir/master_config")
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

  //CLOSE MESSAGE
  println "Results location: ${ outdir }"
}

workflow.onError {
    error_log = file("${master_log_dir}/master_error.log")
    def error_log_msg="""\
    ### ERROR MESSAGE ###
    ${workflow.errorMessage}

    ### ERROR REPORT ###
    ${workflow.errorReport}
    """
    .stripIndent()
    error_log.append(error_log_msg)
}