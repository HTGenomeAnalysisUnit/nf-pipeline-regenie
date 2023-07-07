//Check accessory scripts
regenie_validate_input_java = file("$projectDir/bin/RegenieValidateInput.java", checkIfExists: true)

include { CACHE_JBANG_SCRIPTS       } from '../modules/local/cache_jbang_scripts'
include { VALIDATE_PHENOTYPES       } from '../modules/local/validate_phenotypes'
include { VALIDATE_COVARIATS        } from '../modules/local/validate_covariates'
include { SETUP_MULTIPLE_RUNS       } from '../modules/local/setup_multiple_runs'

workflow PREPARE_PROJECT {
    main:
    if (params.models_table) {
        log.info "Using multi-model mode with models table ${params.models_table}"

        pheno_chunker_r = file("$projectDir/bin/pheno_chunker.R", checkIfExists: true)
        prepare_projects_py = file("$projectDir/bin/prepare_projects.py", checkIfExists: true)
        //conf_template = file("$projectDir/templates/run_parameters.conf", checkIfExists: true)
        //shared_config_file = file(params.shared_config_file, checkIfExists: true)
        models_table = file(params.models_table, checkIfExists: true)
        traits_table = file(params.phenotypes_filename, checkIfExists: true)
        fam_file = file("${params.genotypes_array}.fam", checkIfExists: true)
        
        models_table.copyTo("${params.outdir}/models.tsv")
        
        // This runs phenotyper and prepares the project configuration
        SETUP_MULTIPLE_RUNS(pheno_chunker_r, prepare_projects_py, traits_table, models_table, fam_file) 
        
        phenotype_data = SETUP_MULTIPLE_RUNS.out.analysis_config
            .map{ it.splitCsv(sep: '\t', skip: 1) }
            .toList().transpose()
            .map{ row -> tuple(
                row[0][0],
                file(row[0][1]),
                [
                    cols: row[0][2],
                    binary: {row[0][3] == 'True' ? true : false},
                    model: row[0][4]
                ]
            )} 


        covariate_data = SETUP_MULTIPLE_RUNS.out.analysis_config
            .map{ it.splitCsv(sep: '\t', skip: 1) }
            .toList().transpose()
            .map{ row -> tuple(
                row[0][0],
                file(row[0][5]),
                [
                    cols: row[0][6],
                    cat_cols: row[0][7]
                ]
            )} 

    } else {
        
        //==== READ INPUTS FROM PARAMS ====
        //Phenotypes
        phenotype_data = Channel.of([
                params.project,
                file(params.phenotypes_filename, checkIfExists: true),
                [ 
                    cols: params.phenotypes_columns, 
                    binary: params.phenotypes_binary_trait,
                    model: params.regenie_gwas_test
                ]
            ])

        //Covariates
        if (params.covariates_filename == 'NO_COV_FILE') {
            covar_tmp_file = file("${workflow.workDir}/NO_COV_FILE")
            covar_tmp_file.append('')
            covariates_file = file("$covar_tmp_file")
        } else {
            covariates_file = file(params.covariates_filename, checkIfExists: true)
        }
        covariate_data = Channel.of([
            params.project,
            covariates_file,
            [ 
                cols: params.covariates_columns, 
                cat_cols: params.covariates_cat_columns
            ]
        ])

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