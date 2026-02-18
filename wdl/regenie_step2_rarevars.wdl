version 1.0

## ============================================================================
## Subworkflow: REGENIE Step 2 for Rare Variants
## Handles scatter over chromosomes/gene chunks, then gather results by phenotype
## ============================================================================

import "structs.wdl" as Structs
import "tasks.wdl" as Tasks

workflow RegenieStep2Rarevars {
    input {
        # Project data
        String project_id
        File phenotypes_file
        PhenoMeta pheno_meta
        File covariates_file
        CovarMeta covar_meta
        AccessoryFiles accessory_files
        Array[File] step1_predictions

        # Genotype data for rare vars (one entry per chromosome or single file)
        Array[GenotypeFiles] genotype_files

        # Gene chunk data
        Boolean step2_rarevar_split
        Array[Array[File]]? gene_chunk_files_per_genotype  # gene_chunk_files_per_genotype[i] = gene subset files for genotype_files[i]

        # Rare variant specific files
        File rarevars_set_list
        File rarevars_anno_file
        File rarevars_mask_file

        # Regenie step2 params
        Array[String] chromosomes_list
        String genotypes_rarevar_format
        Int regenie_bsize_step2 = 400
        String regenie_rarevar_min_mac = "1"
        String rarevars_aaf_bins = "0.01,0.05"
        String rarevars_vc_test = "skat,skato,acatv,acato"
        String? rarevars_joint_test
        String? rarevars_vc_maxAAF
        String? regenie_build_mask
        Boolean rarevars_write_mask_snplist = false
        Boolean phenotypes_delete_missings = false
        Boolean regenie_skip_predictions = false
        Boolean regenie_ref_first_step2 = true
        Boolean regenie_firth = true
        Boolean regenie_firth_approx = true
        String regenie_range = ""
        Int maxCatLevels = 10
        String? additional_geno_format
    }

    # Scatter over genotype files; within each, scatter over gene chunks (if splitting)
    scatter (geno_idx in range(length(genotype_files))) {
        GenotypeFiles geno = genotype_files[geno_idx]

        # Determine the number of runs: one per chunk file when splitting, or one run total
        if (step2_rarevar_split && defined(gene_chunk_files_per_genotype)) {
            Array[File] chunks_for_geno = select_first([gene_chunk_files_per_genotype])[geno_idx]
            Int n_chunk_runs = length(chunks_for_geno)
        }
        Int n_runs = select_first([n_chunk_runs, 1])

        scatter (run_idx in range(n_runs)) {
            # Conditionally set the gene chunk file (File? — defined only when splitting)
            if (defined(chunks_for_geno)) {
                File chunk_for_run = select_first([chunks_for_geno])[run_idx]
            }

            call Tasks.RegenieStep2Rarevars as Step2RarevarRun {
                input:
                    project_id = project_id,
                    phenotypes_file = phenotypes_file,
                    pheno_meta = pheno_meta,
                    covariates_file = covariates_file,
                    covar_meta = covar_meta,
                    accessory_files = accessory_files,
                    step1_predictions = step1_predictions,
                    filename = geno.file_prefix,
                    bed_bgen_pgen = geno.primary_file,
                    bim_bgi_pvar = geno.secondary_file,
                    fam_sample_psam = geno.tertiary_file,
                    chromosome = geno.chromosome,
                    gene_chunk_file = chunk_for_run,
                    task_index = run_idx,
                    chromosomes_list = chromosomes_list,
                    rarevars_set_list = rarevars_set_list,
                    rarevars_anno_file = rarevars_anno_file,
                    rarevars_mask_file = rarevars_mask_file,
                    genotypes_rarevar_format = genotypes_rarevar_format,
                    regenie_bsize_step2 = regenie_bsize_step2,
                    regenie_rarevar_min_mac = regenie_rarevar_min_mac,
                    rarevars_aaf_bins = rarevars_aaf_bins,
                    rarevars_vc_test = rarevars_vc_test,
                    rarevars_joint_test = rarevars_joint_test,
                    rarevars_vc_maxAAF = rarevars_vc_maxAAF,
                    regenie_build_mask = regenie_build_mask,
                    rarevars_write_mask_snplist = rarevars_write_mask_snplist,
                    phenotypes_delete_missings = phenotypes_delete_missings,
                    regenie_skip_predictions = regenie_skip_predictions,
                    regenie_ref_first_step2 = regenie_ref_first_step2,
                    regenie_firth = regenie_firth,
                    regenie_firth_approx = regenie_firth_approx,
                    regenie_range = regenie_range,
                    maxCatLevels = maxCatLevels,
                    additional_geno_format = additional_geno_format
            }
        }
    }

    # Flatten all results and logs
    Array[Array[File]] nested_results = flatten(Step2RarevarRun.regenie_results)
    Array[File] all_results = flatten(nested_results)
    Array[File] all_logs = flatten(Step2RarevarRun.step2_log)

    # Parse step2 logs
    call Tasks.RegenieLogParserStep2 as RarevarLogParser {
        input:
            project_id = project_id,
            regenie_step2_logs = all_logs
    }

    # Concatenate results by phenotype
    Array[String] pheno_list = read_lines(write_lines(flatten([split(pheno_meta.cols, ",")])))

    scatter (phenotype in pheno_list) {
        call Tasks.ConcatStep2Results as ConcatRarevarResults {
            input:
                project_id = project_id,
                phenotype = phenotype,
                regenie_gz_files = all_results,
                rarevar_results = true
        }
    }

    output {
        Array[File] merged_results_by_pheno = ConcatRarevarResults.merged_results
        Array[File] merged_results_tbi = ConcatRarevarResults.merged_results_tbi
        Array[String] phenotypes = pheno_list
        File step2_parsed_log = RarevarLogParser.parsed_log
    }
}
