pheno_chunker_r = file("$projectDir/bin/pheno_chunker.R", checkIfExists: true)
prepare_projects_py = file("$projectDir/bin/prepare_projects.py", checkIfExists: true)
conf_template = file("$projectDir/templates/run_parameters.conf", checkIfExists: true)

shared_config_file = file(params.shared_config_file, checkIfExists: true)
models_table = file(params.models_table, checkIfExists: true)
traits_table = file(params.traits_table, checkIfExists: true)

if(!params.master_outdir) {
    exit 1, "Please provide master output directory setting master_outdir parameter"
}
master_outdir = file(params.master_outdir)
master_outdir.mkdirs()

if (params.master_log_dir) {
    master_log_dir = file(params.master_log_dir)
    master_log_dir.mkdirs()
}

for (p in ["genotypes_imputed", "imputed_snplist"]) {
    file(params[p], checkIfExists: true)
}
for (x in ["fam", "bim", "bed"]) {
    file("${params.genotypes_array}.${x}", checkIfExists: true)
}

fam_file = file("${params.genotypes_array}.fam", checkIfExists: true)

include { PREPARE_PROJECTS } from '../modules/local/prepare_projects' addParams(outdir: "$master_outdir")
include { SUBMIT_GWAS_RUN } from '../modules/local/submit_gwas_run' addParams(outdir: "$master_outdir")

workflow SETUP_MULTIPLE_RUNS {
    PREPARE_PROJECTS(pheno_chunker_r, prepare_projects_py, traits_table, models_table, fam_file, conf_template) // This runs phenotyper and prepares the project configuration
    SUBMIT_GWAS_RUN(PREPARE_PROJECTS.out.chunks.flatten(), shared_config_file) //This spawn a nextflow process for each configured run. Probably better 1 or 2 at a time
}

workflow.onComplete {
    models_table.copyTo("$master_outdir")

    //CLOSE MESSAGE
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Results location: ${ master_outdir }"
}