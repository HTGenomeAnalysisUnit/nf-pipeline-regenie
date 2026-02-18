version 1.0

## ============================================================================
## Subworkflow: Run a Single Project Analysis
## ============================================================================
## Executes the full REGENIE pipeline for ONE project:
##   validate inputs → step 1 → step 2 GWAS → step 2 rarevars
##   → process results → generate reports
##
## Called by the main RegenieGWAS workflow, potentially in a scatter
## when running multi-project or multi-model modes.
## ============================================================================

import "structs.wdl" as Structs
import "tasks.wdl" as Tasks
import "regenie_step1.wdl" as Step1WF
import "regenie_step2_gwas.wdl" as Step2GwasWF
import "regenie_step2_rarevars.wdl" as Step2RarevarWF
import "process_gwas_results.wdl" as ProcessGwasWF
import "process_rarevar_results.wdl" as ProcessRarevarWF

workflow RunSingleAnalysis {
    input {
        ## ========== PROJECT-SPECIFIC INPUTS ==========
        String project_id
        File phenotypes_file
        String phenotypes_columns               # comma-separated column names
        Boolean phenotypes_binary_trait
        String regenie_gwas_test = "additive"   # genetic model

        File? covariates_file
        String covariates_columns = ""
        String covariates_cat_columns = ""

        String? interaction_cov
        String? interaction_snp
        File? condition_list
        File? extract_snps_list
        File? extract_genes_list

        ## ========== SHARED: GENOTYPE DATA ==========
        File genotypes_array_bed
        File genotypes_array_bim
        File genotypes_array_fam

        # Imputed data for GWAS
        Array[GenotypeFiles]? genotypes_imputed
        String genotypes_imputed_format = "bgen"

        # Rare variant data
        Array[GenotypeFiles]? genotypes_rarevar
        String genotypes_rarevar_format = "bgen"

        # Pre-computed GWAS chunks from main workflow (shared across projects)
        Array[Array[String]]? gwas_chunks_per_geno

        # LD panel for clumping
        String ld_panel = "NO_LD_FILE"
        Array[File]? ld_panel_beds
        Array[File]? ld_panel_bims
        Array[File]? ld_panel_fams
        Array[String]? ld_panel_chroms

        ## ========== SHARED: SETTINGS ==========
        String genotypes_build
        Array[String] chromosomes
        String project_date = ""
        Boolean phenotypes_delete_missings = false
        Int maxCatLevels = 10

        # Additional genotype data (conditional/interaction analysis)
        File? additional_geno_bgen_pgen_bed
        File? additional_geno_sample_psam_fam
        File? additional_geno_bgi_pvar_bim
        String? additional_geno_format

        # Rare variant accessory files
        File? rarevar_set_list_file
        File? rarevar_anno_file
        File? rarevar_mask_file

        # Gene annotation
        File? genes_bed
        File? genes_ranges

        # Step1 QC settings
        Boolean perform_step1_qc = true
        String qc_maf = "0.01"
        String qc_mac = "100"
        String qc_geno = "0.05"
        String qc_hwe = "1e-15"
        String qc_mind = "0.1"

        # Step1 pruning settings
        Boolean prune_enabled = false
        Float prune_maf = 0.01
        Int prune_window_kbsize = 1000
        Int prune_step_size = 100
        Float prune_r2_threshold = 0.8

        # Regenie step1 settings
        Int regenie_bsize_step1 = 1000
        Boolean regenie_premade_predictions = false
        Array[File]? premade_prediction_files
        Boolean regenie_skip_predictions = false
        Boolean save_step1_predictions = true
        Boolean regenie_force_step1 = false
        Boolean regenie_ref_first_step1 = false
        Boolean step1_use_loocv = false
        Int step1_niter = 30
        Int step1_n_chunks = 100

        # Regenie step2 settings
        Int regenie_bsize_step2 = 400
        Boolean regenie_ref_first_step2 = true
        String regenie_min_imputation_score = "0.00"
        String regenie_gwas_min_mac = "50"
        Boolean regenie_firth = true
        Boolean regenie_firth_approx = true
        String regenie_range = ""

        # GWAS splitting
        Boolean step2_gwas_split = true

        # Rare variant step2 settings
        Boolean step2_rarevar_split = true
        Int step2_rarevar_chunk_size = 200
        String regenie_rarevar_min_mac = "1"
        String rarevars_aaf_bins = "0.01,0.05"
        String rarevars_vc_test = "skat,skato,acatv,acato"
        String rarevars_joint_test = "minp,acat"
        String rarevars_vc_maxAAF = "0.05"
        String regenie_build_mask = "max"
        Boolean rarevars_write_mask_snplist = false

        # Results annotation
        Float annotation_min_log10p = 7.3
        Int annotation_interval_kb = 25
        Boolean clumping = true
        Float clump_p1 = 5e-8
        Float clump_p2 = 1e-4
        Int clump_kb = 250

        # Rare variant results
        Float rarevar_min_log10p = 5.0
        String rarevar_stat_test = "BONF_bygroup"
        Float rarevar_stat_test_threshold = 1.3

        # Reports
        Boolean make_report = true
        String manhattan_annotations = "genes"
        Int regional_plot_window_kb = 300
        Int n_top_loci_plot = 5
        File? gwas_report_template
        File? rarevar_report_template
        File? quarto_report_css

        # Bin scripts
        File? collapse_closegenes_py
        String pipeline_version = "v1.9.4"
    }

    ## ========================================================================
    ## Build structured metadata
    ## ========================================================================

    PhenoMeta pheno_meta = object {
        cols: phenotypes_columns,
        binary: phenotypes_binary_trait,
        model: regenie_gwas_test
    }

    CovarMeta covar_meta = object {
        cols: covariates_columns,
        cat_cols: covariates_cat_columns,
        gxe: interaction_cov,
        gxg: interaction_snp
    }

    AccessoryFiles accessory_files = object {
        condition_list: condition_list,
        additional_bgen_pgen_bed: additional_geno_bgen_pgen_bed,
        additional_sample_psam_fam: additional_geno_sample_psam_fam,
        additional_bgi_pvar_bim: additional_geno_bgi_pvar_bim,
        extract_snps_list: extract_snps_list,
        extract_genes_list: extract_genes_list
    }

    ## ========================================================================
    ## VALIDATE INPUTS
    ## ========================================================================

    call Tasks.ValidatePhenotypes {
        input:
            project_id = project_id,
            phenotypes_file = phenotypes_file,
            pheno_meta = pheno_meta
    }

    if (defined(covariates_file) && covariates_columns != "") {
        call Tasks.ValidateCovariates {
            input:
                project_id = project_id,
                covariates_file = select_first([covariates_file]),
                covar_meta = covar_meta,
                accessory_files = accessory_files
        }
    }

    File validated_pheno = ValidatePhenotypes.validated_phenotypes
    File validated_covar = select_first([
        ValidateCovariates.validated_covariates,
        select_first([covariates_file, phenotypes_file])
    ])

    ## ========================================================================
    ## REGENIE STEP 1
    ## ========================================================================

    call Step1WF.RegenieStep1 {
        input:
            project_id = project_id,
            phenotypes_file = validated_pheno,
            pheno_meta = pheno_meta,
            covariates_file = validated_covar,
            covar_meta = covar_meta,
            accessory_files = accessory_files,
            genotyped_bed = genotypes_array_bed,
            genotyped_bim = genotypes_array_bim,
            genotyped_fam = genotypes_array_fam,
            regenie_skip_predictions = regenie_skip_predictions,
            regenie_premade_predictions = regenie_premade_predictions,
            premade_prediction_files = premade_prediction_files,
            perform_step1_qc = perform_step1_qc,
            qc_maf = qc_maf,
            qc_mac = qc_mac,
            qc_geno = qc_geno,
            qc_hwe = qc_hwe,
            qc_mind = qc_mind,
            prune_enabled = prune_enabled,
            prune_maf = prune_maf,
            prune_window_kbsize = prune_window_kbsize,
            prune_step_size = prune_step_size,
            prune_r2_threshold = prune_r2_threshold,
            regenie_bsize_step1 = regenie_bsize_step1,
            step1_n_chunks = step1_n_chunks,
            step1_niter = step1_niter,
            phenotypes_delete_missings = phenotypes_delete_missings,
            regenie_force_step1 = regenie_force_step1,
            regenie_ref_first_step1 = regenie_ref_first_step1,
            step1_use_loocv = step1_use_loocv,
            maxCatLevels = maxCatLevels,
            additional_geno_format = additional_geno_format,
            project_name = project_id
    }

    ## ========================================================================
    ## REGENIE STEP 2 - GWAS (Common Variants)
    ## ========================================================================

    if (defined(genotypes_imputed)) {
        Array[GenotypeFiles] gwas_geno = select_first([genotypes_imputed])

        # Use pre-computed chunks (shared across projects) or default
        scatter (geno_idx_gwas in range(length(gwas_geno))) {
            Array[String] gwas_chunk_for_geno = if (step2_gwas_split && defined(gwas_chunks_per_geno))
                then select_first([gwas_chunks_per_geno])[geno_idx_gwas]
                else ["SINGLE_CHUNK"]
        }

        call Step2GwasWF.RegenieStep2Gwas {
            input:
                project_id = project_id,
                phenotypes_file = validated_pheno,
                pheno_meta = pheno_meta,
                covariates_file = validated_covar,
                covar_meta = covar_meta,
                accessory_files = accessory_files,
                step1_predictions = RegenieStep1.regenie_step1_out,
                genotype_files = gwas_geno,
                step2_gwas_split = step2_gwas_split,
                chunks_per_genotype = gwas_chunk_for_geno,
                chromosomes_list = chromosomes,
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

        # Process GWAS results
        if (defined(collapse_closegenes_py) && defined(genes_bed)) {
            call ProcessGwasWF.ProcessGwasResults {
                input:
                    project_id = project_id,
                    phenotypes = RegenieStep2Gwas.phenotypes,
                    merged_results_by_pheno = RegenieStep2Gwas.merged_results_by_pheno,
                    genes_bed = select_first([genes_bed]),
                    genes_ranges = genes_ranges,
                    annotation_min_log10p = annotation_min_log10p,
                    annotation_interval_kb = annotation_interval_kb,
                    collapse_closegenes_py = select_first([collapse_closegenes_py]),
                    clumping = clumping,
                    clump_p1 = clump_p1,
                    clump_p2 = clump_p2,
                    clump_kb = clump_kb,
                    genotypes_build = genotypes_build,
                    genotypes_imputed_format = genotypes_imputed_format,
                    chromosomes_list = chromosomes,
                    ld_panel = ld_panel,
                    ld_panel_beds = ld_panel_beds,
                    ld_panel_bims = ld_panel_bims,
                    ld_panel_fams = ld_panel_fams,
                    ld_panel_chroms = ld_panel_chroms,
                    processed_gwas_genotypes = gwas_geno
            }
        }

        # Generate GWAS reports
        if (make_report && defined(gwas_report_template) && defined(quarto_report_css)) {
            scatter (report_idx in range(length(RegenieStep2Gwas.phenotypes))) {
                call Tasks.ReportGwas {
                    input:
                        project_id = project_id,
                        phenotype = RegenieStep2Gwas.phenotypes[report_idx],
                        phenotype_file_validated = validated_pheno,
                        covar_metadata = covar_meta,
                        regenie_merged_results = RegenieStep2Gwas.merged_results_by_pheno[report_idx],
                        annotated_tophits = select_first([ProcessGwasResults.annotated_results])[report_idx],
                        annotated_toploci = select_first([ProcessGwasResults.annotloci])[report_idx],
                        report_template = select_first([gwas_report_template]),
                        quarto_report_css = select_first([quarto_report_css]),
                        project_date = project_date,
                        pipeline_version = pipeline_version,
                        manhattan_annotations = manhattan_annotations,
                        annotation_min_log10p = annotation_min_log10p,
                        n_top_loci_plot = n_top_loci_plot,
                        regional_plot_window_kb = regional_plot_window_kb,
                        genotypes_build = genotypes_build
                }
            }
        }
    }

    ## ========================================================================
    ## REGENIE STEP 2 - RARE VARIANTS
    ## ========================================================================

    if (defined(genotypes_rarevar) && defined(rarevar_set_list_file) && defined(rarevar_anno_file) && defined(rarevar_mask_file)) {
        Array[GenotypeFiles] rarevar_geno = select_first([genotypes_rarevar])

        # Compute gene chunks inline (not shared across projects)
        if (step2_rarevar_split) {
            scatter (rv_geno in rarevar_geno) {
                call Tasks.MakeGenesChunks {
                    input:
                        filename = rv_geno.file_prefix,
                        primary_file = rv_geno.primary_file,
                        secondary_file = rv_geno.secondary_file,
                        tertiary_file = rv_geno.tertiary_file,
                        chromosome = rv_geno.chromosome,
                        set_list_file = select_first([rarevar_set_list_file]),
                        step2_rarevar_chunk_size = step2_rarevar_chunk_size,
                        chromosomes = chromosomes
                }
            }
        }

        call Step2RarevarWF.RegenieStep2Rarevars {
            input:
                project_id = project_id,
                phenotypes_file = validated_pheno,
                pheno_meta = pheno_meta,
                covariates_file = validated_covar,
                covar_meta = covar_meta,
                accessory_files = accessory_files,
                step1_predictions = RegenieStep1.regenie_step1_out,
                genotype_files = rarevar_geno,
                step2_rarevar_split = step2_rarevar_split,
                gene_chunk_files_per_genotype = MakeGenesChunks.chunk_files,
                rarevars_set_list = select_first([rarevar_set_list_file]),
                rarevars_anno_file = select_first([rarevar_anno_file]),
                rarevars_mask_file = select_first([rarevar_mask_file]),
                chromosomes_list = chromosomes,
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

        # Process rare variant results
        call ProcessRarevarWF.ProcessRarevarResults {
            input:
                project_id = project_id,
                phenotypes = RegenieStep2Rarevars.phenotypes,
                merged_results_by_pheno = RegenieStep2Rarevars.merged_results_by_pheno,
                rarevar_min_log10p = rarevar_min_log10p
        }

        # Generate rare variant reports
        if (make_report && defined(rarevar_report_template) && defined(quarto_report_css)) {
            scatter (rv_report_idx in range(length(RegenieStep2Rarevars.phenotypes))) {
                call Tasks.ReportRarevar {
                    input:
                        project_id = project_id,
                        phenotype = RegenieStep2Rarevars.phenotypes[rv_report_idx],
                        phenotype_file_validated = validated_pheno,
                        covar_metadata = covar_meta,
                        regenie_merged_results = RegenieStep2Rarevars.merged_results_by_pheno[rv_report_idx],
                        annotated_tophits = ProcessRarevarResults.filtered_results[rv_report_idx],
                        report_template = select_first([rarevar_report_template]),
                        quarto_report_css = select_first([quarto_report_css]),
                        project_date = project_date,
                        pipeline_version = pipeline_version,
                        genotypes_build = genotypes_build,
                        rarevar_min_log10p = rarevar_min_log10p,
                        rarevar_stat_test_threshold = rarevar_stat_test_threshold,
                        rarevar_stat_test = rarevar_stat_test
                }
            }
        }
    }

    ## ========================================================================
    ## OUTPUTS
    ## ========================================================================

    output {
        # Step 1 outputs
        Array[File] step1_predictions = RegenieStep1.regenie_step1_out
        File? step1_log = RegenieStep1.step1_parsed_log

        # Input validation
        File validated_phenotypes = ValidatePhenotypes.validated_phenotypes
        File phenotype_validation_log = ValidatePhenotypes.validation_log
        File? covariate_validation_log = ValidateCovariates.validation_log

        # GWAS outputs
        Array[File]? gwas_merged_results = RegenieStep2Gwas.merged_results_by_pheno
        Array[File]? gwas_merged_results_tbi = RegenieStep2Gwas.merged_results_tbi
        File? gwas_step2_log = RegenieStep2Gwas.step2_parsed_log
        Array[File]? gwas_filtered_results = ProcessGwasResults.filtered_results
        Array[File]? gwas_annotated_results = ProcessGwasResults.annotated_results
        Array[File?]? gwas_toploci = ProcessGwasResults.toploci
        Array[File?]? gwas_annotloci = ProcessGwasResults.annotloci
        Array[File]? gwas_reports = ReportGwas.report_html

        # Rare variant outputs
        Array[File]? rarevar_merged_results = RegenieStep2Rarevars.merged_results_by_pheno
        Array[File]? rarevar_merged_results_tbi = RegenieStep2Rarevars.merged_results_tbi
        File? rarevar_step2_log = RegenieStep2Rarevars.step2_parsed_log
        Array[File]? rarevar_filtered_results = ProcessRarevarResults.filtered_results
        Array[File]? rarevar_processed_results = ProcessRarevarResults.processed_results
        Array[File]? rarevar_reports = ReportRarevar.report_html
    }
}
