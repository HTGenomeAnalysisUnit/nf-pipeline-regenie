pheno_chunker_r = file("$projectDir/bin/pheno_chunker.R", checkIfExists: true)
prepare_projects_py = file("$projectDir/bin/prepare_projects.py", checkIfExists: true)
//conf_template = file("$projectDir/templates/run_parameters.conf", checkIfExists: true)
shared_config_file = file(params.shared_config_file, checkIfExists: true)
models_table = file(params.models_table, checkIfExists: true)
traits_table = file(params.traits_table, checkIfExists: true)

models_table.copyTo("${params.outdir}/models.tsv")

// if (params.master_log_dir == null) {
//     master_log_dir = file("${master_outdir}/execution_log")
// } else {
//     master_log_dir = file(params.master_log_dir)
// }
// master_log_dir.mkdirs()

// for (p in ["genotypes_imputed", "imputed_snplist"]) {
//     file(params[p], checkIfExists: true)
// }
// for (x in ["fam", "bim", "bed"]) {
//     file("${params.genotypes_array}.${x}", checkIfExists: true)
// }

fam_file = file("${params.genotypes_array}.fam", checkIfExists: true)

include { PREPARE_PROJECT_CONFIG  } from '../modules/local/prepare_projects_config'

workflow SETUP_MULTIPLE_RUNS {
    PREPARE_PROJECTS(pheno_chunker_r, prepare_projects_py, traits_table, models_table, fam_file, conf_template) // This runs phenotyper and prepares the project configuration
    SUBMIT_GWAS_RUN(PREPARE_PROJECTS.out.chunks.flatten(), shared_config_file) //This spawn a nextflow process for each configured run. Probably better 1 or 2 at a time
}