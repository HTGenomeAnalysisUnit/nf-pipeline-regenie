process PROCESS_RAREVAR_RESULTS {
  tag "${project_id}_${phenotype}"
  
  publishDir "${params.outdir}/${project_id}/results/${params.rarevar_results ? 'rarevar' : 'gwas'}", mode: 'copy'

  input:
    tuple val(project_id), val(phenotype), path(regenie_result_gz)

  output:
    tuple val(project_id), val(phenotype), path("${regenie_result_gz.baseName}.correctedP.gz"), emit: results_processed

  script:
  """
  process_regenie_rarevar.py $regenie_result_gz
  """
}