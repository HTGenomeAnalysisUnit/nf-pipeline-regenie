version 1.0

## ============================================================================
## Subworkflow: REGENIE Step 2 for GWAS (common variants)
## Handles scatter over chromosomes/chunks, then gather results by phenotype
## ============================================================================

import "structs.wdl" as Structs
import "tasks.wdl" as Tasks

workflow RegenieStep2Gwas {
    input {
        # Project data
        String project_id
        File phenotypes_file
        PhenoMeta pheno_meta
        File covariates_file
        CovarMeta covar_meta
        AccessoryFiles accessory_files
        Array[File] step1_predictions

        # Genotype data (one entry per chromosome or single file)
        Array[GenotypeFiles] genotype_files

        # Chunk data: parallel arrays with genotype_files
        # If step2_gwas_split is true, chunks_per_geno contains the chunk coordinates
        Boolean step2_gwas_split
        Array[Array[String]] chunks_per_genotype  # chunks_per_genotype[i] = list of chunk coords for genotype_files[i], or ["SINGLE_CHUNK"]

        # Regenie step2 params
        Array[String] chromosomes_list
        String genotypes_imputed_format
        Int regenie_bsize_step2 = 400
        String regenie_gwas_min_mac = "50"
        String regenie_min_imputation_score = "0.00"
        Boolean phenotypes_delete_missings = false
        Boolean regenie_skip_predictions = false
        Boolean regenie_ref_first_step2 = true
        Boolean regenie_firth = true
        Boolean regenie_firth_approx = true
        String regenie_range = ""
        Int maxCatLevels = 10
        String? additional_geno_format
        Boolean save_step2_logs = true
    }

    # Flatten genotype_files x chunks into individual step2 calls
    # We scatter over genotype files first, then over chunks within each
    scatter (geno_idx in range(length(genotype_files))) {
        GenotypeFiles geno = genotype_files[geno_idx]
        Array[String] chunks = chunks_per_genotype[geno_idx]

        scatter (chunk in chunks) {
            call Tasks.RegenieStep2Gwas as Step2GwasRun {
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
                    chunk = chunk,
                    chromosomes_list = chromosomes_list,
                    genotypes_imputed_format = genotypes_imputed_format,
                    regenie_bsize_step2 = regenie_bsize_step2,
                    regenie_gwas_min_mac = regenie_gwas_min_mac,
                    regenie_min_imputation_score = regenie_min_imputation_score,
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
    Array[Array[File]] nested_results = flatten(Step2GwasRun.regenie_results)
    Array[File] all_results = flatten(nested_results)
    Array[File] all_logs = flatten(Step2GwasRun.step2_log)

    # Parse step2 logs
    call Tasks.RegenieLogParserStep2 as GwasLogParser {
        input:
            project_id = project_id,
            regenie_step2_logs = all_logs
    }

    # Concatenate results by phenotype
    # In WDL, we need to know the phenotype list to group results
    Array[String] pheno_list = read_lines(write_lines(flatten([split(pheno_meta.cols, ",")])))

    scatter (phenotype in pheno_list) {
        # Filter results files by phenotype name (files contain phenotype in the name)
        # Since WDL doesn't have dynamic grouping, we pass all results and filter in the concat task
        call Tasks.ConcatStep2Results as ConcatGwasResults {
            input:
                project_id = project_id,
                phenotype = phenotype,
                regenie_gz_files = all_results,
                rarevar_results = false
        }
    }

    output {
        Array[File] merged_results_by_pheno = ConcatGwasResults.merged_results
        Array[File] merged_results_tbi = ConcatGwasResults.merged_results_tbi
        Array[String] phenotypes = pheno_list
        File step2_parsed_log = GwasLogParser.parsed_log
    }
}
