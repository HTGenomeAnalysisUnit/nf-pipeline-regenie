
include { CHECK_PROJECT             } from '../modules/local/check_project'
include { VALIDATE_PHENOTYPES       } from '../modules/local/validate_phenotypes'
include { VALIDATE_COVARIATS        } from '../modules/local/validate_covariates'
include { SETUP_MULTIPLE_RUNS       } from '../modules/local/setup_multiple_runs'

workflow PREPARE_PROJECT {
    main:
    //Make dummy files in case some optional files are not provided
    tmp_files = [:]
    for (x in ['COV', 'CONDITION', 'ADDITIONAL_GENO', 'ADDITIONAL_GENO_SAMPLE']) {
        tmp_filename = "${workflow.workDir}/NO_${x}_FILE"
        tmp_files[x] = file(tmp_filename)
        file(tmp_filename).append('')
    }

    //Set additional genotype files for conditional / interaction analysis
    if (!params.additional_geno_file || params.additional_geno_file == 'NA' || params.additional_geno_file == '') {
        additional_bgen_pgen_bed = file(tmp_files['ADDITIONAL_GENO'], checkIfExists: true)
        additional_sample_psam_fam = file("${tmp_files['ADDITIONAL_GENO']}.sample")
        additional_bgi_pvar_bim = file("${tmp_files['ADDITIONAL_GENO']}.bgi")
    } else {
        if(params.additional_geno_format == 'bgen') {
            additional_bgen_pgen_bed = file("${params.additional_geno_file}.bgen", checkIfExists: true)
            additional_sample_psam_fam = file("${params.additional_geno_file}.sample", checkIfExists: true)
            additional_bgi_pvar_bim = file("${params.additional_geno_file}.bgen.bgi", checkIfExists: true)
        }
        if(params.additional_geno_format == 'bed') {
            additional_bgen_pgen_bed = file("${params.additional_geno_file}.bed", checkIfExists: true)
            additional_sample_psam_fam = file("${params.additional_geno_file}.fam", checkIfExists: true)
            additional_bgi_pvar_bim = file("${params.additional_geno_file}.bim", checkIfExists: true)
        }
        if(params.additional_geno_format == 'pgen') {
            additional_bgen_pgen_bed = file("${params.additional_geno_file}.pgen", checkIfExists: true)
            additional_sample_psam_fam = file("${params.additional_geno_file}.psam", checkIfExists: true)
            additional_bgi_pvar_bim = file("${params.additional_geno_file}.pvar", checkIfExists: true)
        }   
    }

    if (params.models_table) {
        log.info "==> INIT RUN: Using multi-model mode with models table ${params.models_table}"

        pheno_chunker_r = file("$projectDir/bin/pheno_chunker.R", checkIfExists: true)
        prepare_projects_py = file("$projectDir/bin/prepare_projects.py", checkIfExists: true)
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
                    binary: "${row[0][3] == 'True' ? true : false}",
                    model: row[0][4]
                ]
            )} 

        covariate_data = SETUP_MULTIPLE_RUNS.out.analysis_config
            .map{ it.splitCsv(sep: '\t', skip: 1) }
            .toList().transpose()
            .map{ row -> tuple(
                row[0][0],
                file("${!row[0][5] || row[0][5] == 'NA' || row[0][5] == '' || row[0][5] == 'NO_COV_FILE' ? tmp_files['COV'] : row[0][5]}", checkIfExists: true),
                [
                    cols: row[0][6],
                    cat_cols: row[0][7],
                    gxe: params.interaction_cov,
                    gxg: params.interaction_snp
                ],
                [
                    file("${!params.condition_list || params.condition_list == 'NA' || params.condition_list == '' ? tmp_files['CONDITION'] : params.condition_list}", checkIfExists: true),
                    additional_bgen_pgen_bed,
                    additional_sample_psam_fam,
                    additional_bgi_pvar_bim
                ]
            )} 

    } else if (params.projects_table) {
        log.info "==> INIT RUN: Using multi-project mode with projects table ${params.projects_table}"

        //header: ['project_id','pheno_file','pheno_cols','pheno_binary','pheno_model','cov_file','cov_cols','cov_cat_cols']
        phenotype_data = Channel.fromPath(params.projects_table, checkIfExists: true)
            .splitCsv(sep: '\t', header: true)
            .map{ row -> tuple(
                row['project_id'],
                file(row['pheno_file']),
                [
                    cols: row['pheno_cols'],
                    binary: "${row['pheno_binary'] == 'True' ? true : false}",
                    model: row['pheno_model']
                ]
            )} 

        covariate_data = Channel.fromPath(params.projects_table, checkIfExists: true)
            .splitCsv(sep: '\t', header: true)
            .map{ row -> tuple(
                row['project_id'],
                file("${!row['cov_file'] || row['cov_file'] == 'NA' || row['cov_file'] == '' || row['cov_file'] == 'NO_COV_FILE' ? tmp_files['COV'] : row['cov_file']}", checkIfExists: true),
                [
                    cols: row['cov_cols'],
                    cat_cols: row['cov_cat_cols'],
                    gxe: row['interaction_cov'],
                    gxg: row['interaction_snp']
                ],
                [
                    file("${!row['condition_list'] || row['condition_list'] == 'NA' || row['condition_list'] == '' ? tmp_files['CONDITION'] : row['condition_list']}", checkIfExists: true),
                    additional_bgen_pgen_bed,
                    additional_sample_psam_fam,
                    additional_bgi_pvar_bim
                ]
            )} 
    } else {
        log.info "==> INIT RUN: Using single run execution"

        //==== READ INPUTS FROM PARAMS ====
        //Phenotypes
        phenotype_data = Channel.of([
                params.project,
                file(params.phenotypes_filename, checkIfExists: true),
                [ 
                    cols: params.phenotypes_columns, 
                    binary: "${params.phenotypes_binary_trait}",
                    model: params.regenie_gwas_test
                ]
            ])

        //Covariates
        covariate_data = Channel.of([
            params.project,
            file("${!params.covariates_filename || params.covariates_filename == 'NA' || params.covariates_filename == '' || params.covariates_filename == 'NO_COV_FILE' ? tmp_files['COV'] : params.covariates_filename}", checkIfExists: true),
            [ 
                cols: params.covariates_columns, 
                cat_cols: params.covariates_cat_columns,
                gxe: params.interaction_cov,
                gxg: params.interaction_snp
            ],
            [
                file("${!params.condition_list || params.condition_list == 'NA' || params.condition_list == '' ? tmp_files['CONDITION'] : params.condition_list}", checkIfExists: true),
                additional_bgen_pgen_bed,
                additional_sample_psam_fam,
                additional_bgi_pvar_bim
            ]
        ])

    }

    //==== VALIDATE PHENOTYPE INPUT ====
    VALIDATE_PHENOTYPES (
        phenotype_data
    )

    //==== VALIDATE COVARIATE INPUT ====
    covariate_data.branch {
        with_covars: it[1].name != 'NO_COV_FILE'
        no_covars: true
    }.set { validate_covars_inputs }

    VALIDATE_COVARIATS (
        validate_covars_inputs.with_covars
    )

    validated_covars = VALIDATE_COVARIATS.out.covariates_file_validated
        .mix(validate_covars_inputs.no_covars)
    
    validated_covars_logs = VALIDATE_COVARIATS.out.covariates_file_validated_log
        .mix(validate_covars_inputs.no_covars.map { tuple(it[0], file("NO_COV_LOG")) } )

    //==== ASSEMBLE OUTPUT CHANNELS ====
    project_data = VALIDATE_PHENOTYPES.out.phenotypes_file_validated
        .join(validated_covars)
    
    input_validation_logs = VALIDATE_PHENOTYPES.out.phenotypes_file_validated_log
        .join(validated_covars_logs)

    CHECK_PROJECT(project_data.count(), project_data)
    
    emit:
    project_data //[project_id, pheno_file, pheno_meta(cols, binary, model), covar_file, covar_meta(cols, cat_cols), [accessory_files]]
    input_validation_logs //[project_id, pheno_log_file, covar_log_file]
}
