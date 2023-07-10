process REGENIE_LOG_PARSER_STEP2 {
  label 'small_task'
  
  publishDir "${params.outdir}/${project_id}/logs", mode: 'copy'

  input:
    tuple val(project_id), path(regenie_step2_logs)
    path regenie_log_parser_jar

  output:
    tuple val(project_id), path("${project_id}.step2.log"), emit: regenie_step2_parsed_logs

  """
  java -jar ${regenie_log_parser_jar} ${regenie_step2_logs} --output ${project_id}.step2.log
  """

}
