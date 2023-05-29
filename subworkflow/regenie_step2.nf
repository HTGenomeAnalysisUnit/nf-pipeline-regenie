include { REGENIE_LOG_PARSER_STEP2 as GWAS_LOG_PARSER_STEP2     } from '../modules/local/regenie_log_parser_step2'  addParams(outdir: "${params.outdir}/logs")
include { REGENIE_LOG_PARSER_STEP2 as RAREVAR_LOG_PARSER_STEP2  } from '../modules/local/regenie_log_parser_step2'  addParams(outdir: "${params.outdir}/logs")
include { CONCAT_STEP2_RESULTS as CONCAT_GWAS_RESULTS           } from '../modules/local/concat_step2_results'      addParams(outdir: "${params.outdir}/results", rarevar_results: false)
include { CONCAT_STEP2_RESULTS as CONCAT_RAREVAR_RESULTS        } from '../modules/local/concat_step2_results'      addParams(outdir: "${params.outdir}/results", rarevar_results: true)
include { REGENIE_STEP2_GWAS                                    } from '../modules/local/regenie_step2'             addParams(outdir: "${params.outdir}/logs", save_step2_logs: params.save_step2_logs, regenie_ref_first: params.regenie_ref_first_step2)
include { REGENIE_STEP2_RAREVARS                                } from '../modules/local/regenie_step2'             addParams(outdir: "${params.outdir}/logs", save_step2_logs: params.save_step2_logs, regenie_ref_first: params.regenie_ref_first_step2)

workflow REGENIE_STEP2_WF {
    take:
        processed_genotypes_ch //val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom), val(chunk)
        regenie_step1_out_ch
        phenotypes_file_validated
        covariates_file_validated
        regenie_log_parser_jar

    main:
    //==== INITIALIZATION ====
    step2_log = Channel.empty()
    step2_results_by_pheno = Channel.empty()

    //==== REGENIE STEP 2 FOR COMMON VARS ====
    if (params.run_gwas) {
        REGENIE_STEP2_GWAS (
            regenie_step1_out_ch.collect(),
            processed_genotypes_ch,
            phenotypes_file_validated,
            covariates_file_validated
        )

        //Parse log file
        GWAS_LOG_PARSER_STEP2 (
            REGENIE_STEP2_GWAS.out.regenie_step2_out_log.first(),
            regenie_log_parser_jar
        )

        //Concatenate results
        concat_gwas_results_ch = REGENIE_STEP2_GWAS.out.regenie_step2_out.transpose()
            .map{ tuple(it[3].baseName.replace("${it[0]}_${it[1]}_${it[2]}_",'').replace('.regenie',''), it[3]) }
            .groupTuple().map{ tuple(it[0], it[1].flatten()) }
        
        CONCAT_GWAS_RESULTS(concat_gwas_results_ch)

        step2_log = GWAS_LOG_PARSER_STEP2.out.regenie_step2_parsed_logs
        step2_results_by_pheno = CONCAT_GWAS_RESULTS.out.regenie_results_gz
    }

    //==== REGENIE STEP 2 FOR RARE VARS ====
    if (params.run_rarevar) {
        rarevars_set_list = file(params.rarevar_set_list_file, checkIfExists: true)
        rarevars_anno_file = file(params.rarevar_anno_file, checkIfExists: true)
        rarevars_mask_file = file(params.rarevar_mask_file, checkIfExists: true)

        REGENIE_STEP2_RAREVARS (
            regenie_step1_out_ch.collect(),
            processed_genotypes_ch,
            phenotypes_file_validated,
            covariates_file_validated,
            rarevars_anno_file,
            rarevars_set_list,
            rarevars_mask_file
        )

        //Parse log file
        RAREVAR_LOG_PARSER_STEP2 (
            REGENIE_STEP2_RAREVARS.out.regenie_step2_out_log.first(),
            regenie_log_parser_jar
        )
        //Concatenate results
        concat_rarevar_results_ch = REGENIE_STEP2_RAREVARS.out.regenie_step2_out.transpose()
            .map{ tuple(it[3].baseName.replace("${it[0]}_${it[1]}_${it[2]}_",'').replace('.regenie',''), it[3]) }
            .groupTuple().map{ tuple(it[0], it[1].flatten()) }
        
        CONCAT_RAREVAR_RESULTS(concat_rarevar_results_ch)

        step2_log = RAREVAR_LOG_PARSER_STEP2.out.regenie_step2_parsed_logs
        step2_results_by_pheno = CONCAT_RAREVAR_RESULTS.out.regenie_results_gz
    }
    
    emit:
        regenie_results = step2_results_by_pheno
        regenie_log = step2_log
}