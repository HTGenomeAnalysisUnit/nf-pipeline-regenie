include { QC_FILTER_GENOTYPED                  } from '../modules/local/qc_filter_genotyped'        addParams(logdir: "${params.outdir}/logs")
include { PRUNE_GENOTYPED                      } from '../modules/local/prune_genotyped'            addParams(logdir: "${params.outdir}/logs")
include { REGENIE_STEP1_SPLIT as REGENIE_STEP1 } from '../modules/local/regenie_step1'              addParams(outdir: "${params.outdir}/regenie_step1_preds", logdir: "${params.outdir}/logs", save_step1_predictions: params.save_step1_predictions, use_loocv: params.step1_use_loocv, niter: params.step1_niter, regenie_ref_first: params.regenie_ref_first_step1)
include { REGENIE_LOG_PARSER_STEP1             } from '../modules/local/regenie_log_parser_step1'   addParams(logdir: "${params.outdir}/logs")

workflow REGENIE_STEP1_WF {
    take:
        genotyped_plink_ch
        phenotypes_file_validated
        covariates_file_validated
        regenie_log_parser_jar
  
    main:
    //==== PREPARE GENOTYPE DATA FOR STEP1 ====
    QC_FILTER_GENOTYPED (
        genotyped_plink_ch, phenotypes_file_validated
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
        regenie_step1_out_ch = Channel.fromPath("NO_PREDICTIONS")
        regenie_step1_parsed_logs_ch = Channel.fromPath("NO_LOG")
    } else if (params.regenie_premade_predictions) {
        /* 
        You can load pre-made regenie level 1 preds. 
        You must specify a path of /my/path/regenie_step1_out* 
        A file named regenie_step1_out_pred.list must be present together with files like
        regenie_step1_out_1.loco.gz, regenie_step1_out_2.loco.gz, ... (one per phenotype)
        These can be used in STEP 2 only if phenotypes and covariates are exactly the same 
        used to generate step1 predictions (both files and column designation must match exactly)
        */
        regenie_step1_out_ch = Channel
            .fromPath(params.regenie_premade_predictions, checkIfExists: true)
        regenie_step1_parsed_logs_ch = Channel.fromPath("NO_LOG")
    } else {
        REGENIE_STEP1 (
            genotyped_final_ch,
            QC_FILTER_GENOTYPED.out.genotyped_filtered_snplist_ch,
            QC_FILTER_GENOTYPED.out.genotyped_filtered_id_ch,
            phenotypes_file_validated,
            covariates_file_validated
        )

        REGENIE_LOG_PARSER_STEP1 (
            REGENIE_STEP1.out.regenie_step1_out_log,
            regenie_log_parser_jar
        )

        regenie_step1_out_ch = REGENIE_STEP1.out.regenie_step1_out
        regenie_step1_parsed_logs_ch = REGENIE_LOG_PARSER_STEP1.out.regenie_step1_parsed_logs
    }

    emit:
        regenie_step1_out = regenie_step1_out_ch
        regenie_step1_parsed_logs = regenie_step1_parsed_logs_ch
}