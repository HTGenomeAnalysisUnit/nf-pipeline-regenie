process OPENING_LOG {
	input:
	val(project_data) // [project_id, pheno_file, pheno_meta(cols, binary, model), covar_file, covar_meta(cols, cat_cols)]

	exec:
log.info """
==========================================================
  PROJECT ${project_data[0]} - NF PIPELINE    
==========================================================
PROJECT ID      	: ${project_data[0]}
OUTDIR          	: ${params.outdir}/${project_data[0]}
PHENOTYPE FILE  	: ${project_data[1]}
PHENOTYPES      	: ${project_data[2].cols}
BINARY PHENOTYPES	: ${project_data[2].binary}
MODEL  			 	: ${project_data[2].model}
COVAR FILE         	: ${project_data[3]}
COVARIATES         	: ${project_data[4].cols}
CATEGORICAL COVARS 	: ${project_data[4].cat_cols}

==========================================================
"""
}