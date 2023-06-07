process VALIDATE_PHENOTYPES {
  label 'small_task'
  
  publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*log'
  publishDir "${params.outdir}/validated_input/", mode: 'copy', pattern: '*validated.txt'

  input:
    path phenotypes_file
    path regenie_validate_input_jar

  output:
    path "${phenotypes_file.baseName}.pheno.validated.txt", emit: phenotypes_file_validated
    path "${phenotypes_file.baseName}.pheno.validated.log", emit: phenotypes_file_validated_log

  """
  java -jar ${regenie_validate_input_jar}  --input ${phenotypes_file} --output  ${phenotypes_file.baseName}.pheno.validated.txt --type phenotype
  """
  }
