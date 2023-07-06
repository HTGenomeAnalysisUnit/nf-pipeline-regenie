process REGENIE_LOG_PARSER_STEP1 {
  label 'small_task'
  
  publishDir "${params.logdir}", mode: 'copy'

  input:
    tuple val(project_id), path(regenie_step1_log)
    path regenie_log_parser_jar

  output:
    tuple val(project_id), path("${params.project}.step1.log"), emit: regenie_step1_parsed_logs

  """
  java -jar ${regenie_log_parser_jar} ${regenie_step1_log} --output ${params.project}.step1.log
  """
  }
