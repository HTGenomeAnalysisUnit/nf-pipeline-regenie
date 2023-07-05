process CACHE_JBANG_SCRIPTS {
  label 'small_task'
  
  input:
    path java_script

  output:
    path "compiled.jar", emit: compiled_jar

  """
  mkdir -p jbang_cache
  export JBANG_CACHE_DIR=jbang_cache
  jbang export portable -O=compiled.jar ${java_script}
  """

}
