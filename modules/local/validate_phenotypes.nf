process VALIDATE_PHENOTYPES {
  label 'small_task'
  
  publishDir {"${params.outdir}/${project_id}/logs"}, mode: 'copy', pattern: '*log'
  publishDir {"${params.outdir}/${project_id}/validated_input"}, mode: 'copy', pattern: '*validated.txt'

  input:
    tuple val(project_id), file(phenotypes_file), val(pheno_meta)

  output:
    tuple val(project_id), file("${phenotypes_file.baseName}.pheno.validated.txt"), val(pheno_meta), emit: phenotypes_file_validated
    tuple val(project_id), path("${phenotypes_file.baseName}.pheno.validated.log"), emit: phenotypes_file_validated_log

  """
  RegenieValidateInput.py --input ${phenotypes_file} --output ${phenotypes_file.baseName}.pheno.validated.txt --type phenotype
  """
  }
