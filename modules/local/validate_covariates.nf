process VALIDATE_COVARIATS {
  label 'small_task'
  
  publishDir {"${params.outdir}/${project_id}/logs"}, mode: 'copy', pattern: '*log'
  publishDir {"${params.outdir}/${project_id}/validated_input"}, mode: 'copy', pattern: '*validated.txt'

  input:
    tuple val(project_id), file(covariates_file), val(covar_meta), file(accessory_files)

  output:
    tuple val(project_id), path("${covariates_file.baseName}.cov.validated.txt"), val(covar_meta), file(accessory_files), emit: covariates_file_validated
    tuple val(project_id), path("${covariates_file.baseName}.cov.validated.log"), emit: covariates_file_validated_log

  """
  RegenieValidateInput.py --input ${covariates_file} --output  ${covariates_file.baseName}.cov.validated.txt --type covariate
  """
  }
