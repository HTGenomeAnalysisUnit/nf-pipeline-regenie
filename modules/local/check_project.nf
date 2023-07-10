process CHECK_PROJECT {
	input:
	val(project_size)
	tuple val(project_id), path(pheno_file), val(pheno_meta), path(covar_file), val(covar_meta)

	exec:
	//Check specified regenie test is allowed
	def allowed_tests = ['additive', 'dominant', 'recessive']
  	if (!(pheno_meta.model in allowed_tests)){
    	log.error "Test ${pheno_meta.model} not supported. Allowed tests are $allowed_tests."
    	exit 1
  	}

	if (!pheno_file.exists()) {
		log.error "Phenotype file ${pheno_file} does not exist."
		exit 1
	}

	if (!covar_file.exists()) {
		log.error "Covariate file ${covar_file} does not exist."
		exit 1
	}

log.info """
Successfully checked project ${project_id} configuration.
"""
}