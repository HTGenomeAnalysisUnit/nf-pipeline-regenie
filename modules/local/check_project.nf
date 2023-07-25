process CHECK_PROJECT {
	executor 'local'
	
	input:
	val(project_size)
	tuple val(project_id), path(pheno_file), val(pheno_meta), path(covar_file), val(covar_meta)

	exec:
	pipeline_log_dir = file("${params.logdir}/analysis_config")
  	pipeline_log_dir.mkdirs()

	def msg="""\
	> Project: ${project_id}
		- Phenotype file: ${pheno_file}
		- Phenotype columns: ${pheno_meta.cols}
		- Phenotype model: ${pheno_meta.model}
		- Covariate file: ${covar_file}
		- Covariate columns: ${covar_meta.cols}

	"""
	.stripIndent()

	//Check specified regenie test is allowed
	def allowed_tests = ['additive', 'dominant', 'recessive']
  	if (!(pheno_meta.model in allowed_tests)){
    	log.error "Test ${pheno_meta.model} not supported in project ${project_id}. Allowed tests are $allowed_tests."
    	exit 1
  	}

	log.info "Project ${project_id} initialized."

	projects_config_log = file("${pipeline_log_dir}/projects_configuration.log")
  	projects_config_log.append(msg)
}