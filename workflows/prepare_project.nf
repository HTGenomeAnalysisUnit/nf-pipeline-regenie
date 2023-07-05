//Check accessory scripts
regenie_validate_input_java = file("$projectDir/bin/RegenieValidateInput.java", checkIfExists: true)

include { CACHE_JBANG_SCRIPTS         } from '../modules/local/cache_jbang_scripts'
include { VALIDATE_PHENOTYPES         } from '../modules/local/validate_phenotypes' addParams(outdir: "$outdir")
include { VALIDATE_COVARIATS          } from '../modules/local/validate_covariates' addParams(outdir: "$outdir")

workflow PREPARE_PROJECT {
    main:
    if (params.with_master) {
        log.info "Using multi-model mode with a master table"
        SETUP_MULTIPLE_RUNS()
    } else {
        
        //==== READ INPUTS FROM PARAMS ====
        //Phenotypes
        def phenotype_data = [
                params.project
                file(params.phenotypes_filename, checkIfExists: true),
                [ 
                    cols: params.phenotypes_columns, 
                    binary: params.phenotypes_binary_trait
                    model: params.regenie_gwas_test
                ]
            ]

        //Covariates
        if (params.covariates_filename == 'NO_COV_FILE') {
            covar_tmp_file = file("${workflow.workDir}/NO_COV_FILE")
            covar_tmp_file.append('')
            covariates_file = file("$covar_tmp_file")
        } else {
            covariates_file = file(params.covariates_filename, checkIfExists: true)
        }
        def covariate_data = [
                params.project
                covariates_file,
                [ 
                cols: params.covariates_columns, 
                cat_cols: params.covariates_cat_columns
                ]
            ]

    }

    //==== PREPARE SCRIPTS ====
    CACHE_JBANG_SCRIPTS (
        regenie_validate_input_java
    )

    //==== VALIDATE PHENOTYPE INPUT ====
    VALIDATE_PHENOTYPES (
        phenotype_data,
        CACHE_JBANG_SCRIPTS.out.compiled_jar
    )

    //==== VALIDATE COVARIATE INPUT ====
    covariate_data.branch {
        with_covars: it[3] != 'NO_COV_FILE'
        no_covars: true
    }.set { validate_covars_inputs }

    VALIDATE_COVARIATS (
        validate_covars_inputs.with_covars,
        CACHE_JBANG_SCRIPTS.out.compiled_jar
    )

    validated_covars = VALIDATE_COVARIATS.out.covariates_file_validated
        .mix(validate_covars_inputs.no_covars)
    
    validated_covars_logs = VALIDATE_COVARIATS.out.covariates_file_validated_log
        .mix(validate_covars_inputs.no_covars.map { tuple(it[0], path("NO_COV_LOG")) } )

    
    //==== ASSEMBLE OUTPUT CHANNELS ====
    project_data = VALIDATE_PHENOTYPES.out.phenotypes_file_validated
        .join(validated_covars)
    
    input_validation_logs = VALIDATE_PHENOTYPES.out.phenotypes_file_validated_log
        .join(validated_covars_logs)

    project_data.view()
    input_validation_logs.view()

    emit:
    project_data //[project_id, pheno_file, pheno_meta(cols, binary, model), covar_file, covar_meta(cols, cat_cols)]
    input_validation_logs //[project_id, pheno_log_file, covar_log_file]
}