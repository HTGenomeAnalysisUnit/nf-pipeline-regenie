process REGENIE_LOG_PARSER_STEP2 {
  label 'small_task'
  
  publishDir "${params.logdir}/${project_id}/logs", mode: 'copy'

  input:
    tuple val(project_id), path(regenie_step2_logs)

  output:
    tuple val(project_id), path("${project_id}.step2.log"), emit: regenie_step2_parsed_logs

  """
  RegenieLogParser.py ${regenie_step2_logs[0]} --output ${project_id}.step2.log
  """

}
