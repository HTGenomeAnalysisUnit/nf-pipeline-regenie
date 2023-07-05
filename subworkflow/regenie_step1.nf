regenie_log_parser_java  = file("$projectDir/bin/RegenieLogParser.java", checkIfExists: true)

include { QC_FILTER_GENOTYPED                  } from '../modules/local/qc_filter_genotyped'        addParams(logdir: "${params.outdir}")
include { PRUNE_GENOTYPED                      } from '../modules/local/prune_genotyped'            addParams(logdir: "${params.outdir}")
include { REGENIE_STEP1_SPLIT as REGENIE_STEP1 } from '../modules/local/regenie_step1'              addParams(outdir: "${params.outdir}", logdir: "${params.outdir}", save_step1_predictions: params.save_step1_predictions, use_loocv: params.step1_use_loocv, niter: params.step1_niter, regenie_ref_first: params.regenie_ref_first_step1)
include { REGENIE_LOG_PARSER_STEP1             } from '../modules/local/regenie_log_parser_step1'   addParams(logdir: "${params.outdir}")

workflow REGENIE_STEP1_WF {
    take:
        genotyped_plink_ch //[project_id, bed_file, bim_file, fam_file]
        project_data //[project_id, pheno_file, pheno_meta(cols, binary, model), covar_file, covar_meta(cols, cat_cols)]
  
    main:
    CACHE_JBANG_SCRIPTS (
        regenie_log_parser_java
    )

    //==== PREPARE GENOTYPE DATA FOR STEP1 ====
    qc_input_ch = genotyped_plink_ch.
        .join(project_data.map { tuple(it[0], it[1]) })
    
    QC_FILTER_GENOTYPED (
        qc_input_ch
    )

    if(params.prune_enabled) {
        PRUNE_GENOTYPED (
            QC_FILTER_GENOTYPED.out.genotyped_filtered_files_ch
        )

        genotyped_final_ch = PRUNE_GENOTYPED.out.genotypes_pruned_ch
    } else {
        //no pruning applied, set QCed directly to genotyped_final_ch
        genotyped_final_ch = QC_FILTER_GENOTYPED.out.genotyped_filtered_files_ch
    }

    //==== REGENIE STEP 1 ====
    if (params.regenie_skip_predictions) {
        //skip step 1 predictions completely
        regenie_step1_out_ch = project_data.map{ it + path("NO_PREDICTIONS") }
        regenie_step1_parsed_logs_ch = project_data.map{ tuple(it[0], path("NO_LOG")) }
    } else if (params.regenie_premade_predictions) {
        /* 
        You can load pre-made regenie level 1 preds. 
        You must specify a path of /my/path/regenie_step1_out* 
        A file named regenie_step1_out_pred.list must be present together with files like
        regenie_step1_out_1.loco.gz, regenie_step1_out_2.loco.gz, ... (one per phenotype)
        These can be used in STEP 2 only if phenotypes and covariates are exactly the same 
        used to generate step1 predictions (both files and column designation must match exactly)
        */
        regenie_step1_out_ch =  project_data.map{ it + path(params.regenie_premade_predictions, checkIfExists: true) }
        regenie_step1_parsed_logs_ch = project_data.map{ tuple(it[0], path("NO_LOG")) }
    } else {

        step1_input_ch = project_data.join(genotyped_final_ch)

        REGENIE_STEP1 (step1_input_ch)

        REGENIE_LOG_PARSER_STEP1 (
            REGENIE_STEP1.out.regenie_step1_out_log,
            CACHE_JBANG_SCRIPTS.out.compiled_jar
        )

        regenie_step1_out_ch = REGENIE_STEP1.out.regenie_step1_out
        regenie_step1_parsed_logs_ch = REGENIE_LOG_PARSER_STEP1.out.regenie_step1_parsed_logs
    }

    emit:
        regenie_step1_out = regenie_step1_out_ch
        regenie_step1_parsed_logs = regenie_step1_parsed_logs_ch
}