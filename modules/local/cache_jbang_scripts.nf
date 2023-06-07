process CACHE_JBANG_SCRIPTS {
  label 'small_task'
  
  input:
    path regenie_log_parser_java
    path regenie_filter_java
    path regenie_validate_input_java

  output:
    path "RegenieLogParser.jar", emit: regenie_log_parser_jar
    path "RegenieFilter.jar", emit: regenie_filter_jar
    path "RegenieValidateInput.jar", emit: regenie_validate_input_jar

  """
  mkdir -p jbang_cache
  export JBANG_CACHE_DIR=jbang_cache
  jbang export portable -O=RegenieLogParser.jar ${regenie_log_parser_java}
  jbang export portable -O=RegenieFilter.jar ${regenie_filter_java}
  jbang export portable -O=RegenieValidateInput.jar ${regenie_validate_input_java}
  """

}
