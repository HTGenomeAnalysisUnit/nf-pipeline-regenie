version 1.0

## ============================================================================
## Subworkflow: REGENIE Step 1 (QC, optional pruning, split-L0/L1 regenie)
## ============================================================================

import "structs.wdl" as Structs
import "tasks.wdl" as Tasks

workflow RegenieStep1 {
    input {
        # Project data
        String project_id
        File phenotypes_file
        PhenoMeta pheno_meta
        File covariates_file
        CovarMeta covar_meta
        AccessoryFiles accessory_files

        # Genotyped plink files
        File genotyped_bed
        File genotyped_bim
        File genotyped_fam

        # Step1 control params
        Boolean regenie_skip_predictions = false
        Boolean regenie_premade_predictions = false
        Array[File]? premade_prediction_files

        # QC params
        Boolean perform_step1_qc = true
        String qc_maf = "0.01"
        String qc_mac = "100"
        String qc_geno = "0.05"
        String qc_hwe = "1e-15"
        String qc_mind = "0.1"

        # Pruning params
        Boolean prune_enabled = false
        Float prune_maf = 0.01
        Int prune_window_kbsize = 1000
        Int prune_step_size = 100
        Float prune_r2_threshold = 0.8

        # Regenie step1 params
        Int regenie_bsize_step1 = 1000
        Int step1_n_chunks = 100
        Int step1_niter = 30
        Boolean phenotypes_delete_missings = false
        Boolean regenie_force_step1 = false
        Boolean regenie_ref_first_step1 = false
        Boolean step1_use_loocv = false
        Int maxCatLevels = 10
        String? additional_geno_format
        String project_name
    }

    # If skip predictions, return empty
    if (!regenie_skip_predictions && !regenie_premade_predictions) {

        # === QC Filter ===
        if (perform_step1_qc) {
            call Tasks.QcFilterGenotyped {
                input:
                    project_id = project_id,
                    bed_file = genotyped_bed,
                    bim_file = genotyped_bim,
                    fam_file = genotyped_fam,
                    phenos_tsv = phenotypes_file,
                    qc_maf = qc_maf,
                    qc_mac = qc_mac,
                    qc_geno = qc_geno,
                    qc_hwe = qc_hwe,
                    qc_mind = qc_mind
            }

            # === Optional Pruning ===
            if (prune_enabled) {
                call Tasks.PruneGenotyped {
                    input:
                        project_id = project_id,
                        bed_file = QcFilterGenotyped.qc_bed,
                        bim_file = QcFilterGenotyped.qc_bim,
                        fam_file = QcFilterGenotyped.qc_fam,
                        prune_maf = prune_maf,
                        prune_window_kbsize = prune_window_kbsize,
                        prune_step_size = prune_step_size,
                        prune_r2_threshold = prune_r2_threshold
                }
            }
        }

        # Determine which bed/bim/fam to use for step1
        File step1_bed = select_first([PruneGenotyped.pruned_bed, QcFilterGenotyped.qc_bed, genotyped_bed])
        File step1_bim = select_first([PruneGenotyped.pruned_bim, QcFilterGenotyped.qc_bim, genotyped_bim])
        File step1_fam = select_first([PruneGenotyped.pruned_fam, QcFilterGenotyped.qc_fam, genotyped_fam])

        # === Split L0 ===
        call Tasks.RegenieSplitL0 {
            input:
                project_id = project_id,
                phenotypes_file = phenotypes_file,
                pheno_meta = pheno_meta,
                covariates_file = covariates_file,
                covar_meta = covar_meta,
                accessory_files = accessory_files,
                bed_file = step1_bed,
                bim_file = step1_bim,
                fam_file = step1_fam,
                regenie_bsize_step1 = regenie_bsize_step1,
                step1_n_chunks = step1_n_chunks,
                phenotypes_delete_missings = phenotypes_delete_missings,
                regenie_ref_first = regenie_ref_first_step1,
                maxCatLevels = maxCatLevels,
                additional_geno_format = additional_geno_format
        }

        # === Run L0 jobs in parallel (scatter) ===
        Array[Int] l0_jobs = range(step1_n_chunks)
        scatter (job_idx in l0_jobs) {
            Int job_number = job_idx + 1

            call Tasks.RegenieRunL0 {
                input:
                    job_number = job_number,
                    project_id = project_id,
                    phenotypes_file = phenotypes_file,
                    pheno_meta = pheno_meta,
                    covariates_file = covariates_file,
                    covar_meta = covar_meta,
                    accessory_files = accessory_files,
                    master_file = RegenieSplitL0.master_file,
                    snplists = RegenieSplitL0.snplists,
                    bed_file = step1_bed,
                    bim_file = step1_bim,
                    fam_file = step1_fam,
                    regenie_bsize_step1 = regenie_bsize_step1,
                    phenotypes_delete_missings = phenotypes_delete_missings,
                    regenie_force_step1 = regenie_force_step1,
                    regenie_ref_first = regenie_ref_first_step1,
                    use_loocv = step1_use_loocv,
                    maxCatLevels = maxCatLevels,
                    additional_geno_format = additional_geno_format
            }
        }

        # Flatten all L0 outputs
        Array[File] all_l0_files = flatten(RegenieRunL0.l0_output)

        # === Run L1 jobs per phenotype (scatter) ===
        Array[String] pheno_list = read_lines(write_lines(flatten([split(pheno_meta.cols, ",")])))

        scatter (single_pheno in pheno_list) {
            call Tasks.RegenieRunL1 {
                input:
                    project_id = project_id,
                    single_pheno = single_pheno,
                    phenotypes_file = phenotypes_file,
                    pheno_meta = pheno_meta,
                    covariates_file = covariates_file,
                    covar_meta = covar_meta,
                    accessory_files = accessory_files,
                    master_file = RegenieSplitL0.master_file,
                    snplists = RegenieSplitL0.snplists,
                    runl0_files = all_l0_files,
                    bed_file = step1_bed,
                    bim_file = step1_bim,
                    fam_file = step1_fam,
                    regenie_bsize_step1 = regenie_bsize_step1,
                    niter = step1_niter,
                    phenotypes_delete_missings = phenotypes_delete_missings,
                    regenie_force_step1 = regenie_force_step1,
                    regenie_ref_first = regenie_ref_first_step1,
                    use_loocv = step1_use_loocv,
                    maxCatLevels = maxCatLevels,
                    additional_geno_format = additional_geno_format
            }
        }

        # Concat all pred lists
        call Tasks.ConcatPredLists {
            input:
                project_id = project_id,
                pred_list_files = select_all(RegenieRunL1.pred_list)
        }

        # Parse step1 log (use last pheno's log as representative)
        call Tasks.RegenieLogParserStep1 {
            input:
                project_id = project_id,
                regenie_step1_log = select_first(RegenieRunL1.step1_log),
                project_name = project_name
        }

        # Collect all step1 output files: gz files + merged pred list
        Array[File] all_step1_gz = flatten(select_all(RegenieRunL1.step1_gz_files))
        Array[File] step1_predictions_computed = flatten([all_step1_gz, [ConcatPredLists.merged_pred_list]])
    }

    # Determine final step1 predictions
    Array[File] final_step1_predictions = select_first([
        step1_predictions_computed,
        premade_prediction_files,
        []
    ])

    output {
        Array[File] regenie_step1_out = final_step1_predictions
        File? step1_parsed_log = RegenieLogParserStep1.parsed_log
    }
}
