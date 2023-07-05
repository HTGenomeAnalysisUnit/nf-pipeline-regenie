process VALIDATE_COVARIATS {
  label 'small_task'
  
  publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*log'
  publishDir "${params.outdir}/validated_input/", mode: 'copy', pattern: '*validated.txt'

  input:
    tuple val(project_id), file(covariates_file), val(covar_meta)
    path regenie_validate_input_jar

  output:
    tuple val(project_id), path("${covariates_file.baseName}.cov.validated.txt"), val(covar_meta), emit: covariates_file_validated
    tuple val(project_id), path("${covariates_file.baseName}.cov.validated.log"), emit: covariates_file_validated_log

  """
  java -jar ${regenie_validate_input_jar} --input ${covariates_file} --output  ${covariates_file.baseName}.cov.validated.txt --type covariate
  """
  }
