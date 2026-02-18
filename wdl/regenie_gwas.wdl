version 1.0

## ============================================================================
## Main Workflow: REGENIE GWAS Pipeline (WDL)
## ============================================================================
## Original Nextflow pipeline: nf-pipeline-regenie v1.9.4
## Author: Edoardo Giacopuzzi
## ============================================================================
##
## EXECUTION MODES:
## 1. Single project: Set project, phenotypes_filename, phenotypes_columns, etc.
## 2. Models table:   Set models_table + phenotypes_filename (combined pheno+covar file)
## 3. Projects array: Set projects (Array[ProjectConfig])
##
## The mode is determined automatically:
##   - If models_table is defined → models table mode
##   - Else if projects is defined → projects array mode
##   - Else → single project mode
## ============================================================================

import "structs.wdl" as Structs
import "tasks.wdl" as Tasks
import "run_single_analysis.wdl" as RunAnalysisWF

workflow RegenieGWAS {
    input {
        ## ====================================================================
        ## MODE 1: SINGLE PROJECT INPUTS
        ## ====================================================================
        ## Required only in single project mode. Ignored if models_table or
        ## projects is set.
        String? project                          # Project name
        File? phenotypes_filename                # Phenotype file
        String? phenotypes_columns               # Comma-separated phenotype column names
        Boolean? phenotypes_binary_trait          # true for binary traits

        File? covariates_filename                # Covariate file (optional)
        String covariates_columns = ""           # Comma-separated covariate columns
        String covariates_cat_columns = ""       # Comma-separated categorical covariates

        String? interaction_cov                  # GxE interaction covariate
        String? interaction_snp                  # GxG interaction variant ID
        File? condition_list                     # Variants to condition on
        File? extract_snps_list                  # Restrict GWAS to these variants
        File? extract_genes_list                 # Restrict rare-var analysis to these genes

        ## ====================================================================
        ## MODE 2: MODELS TABLE INPUTS
        ## ====================================================================
        ## When models_table is set, the pipeline runs the phenotype chunker
        ## (R script) to create per-model phenotype/covariate files, then
        ## scatters over models. In this mode, phenotypes_filename must be
        ## the combined file containing all phenotypes AND covariates.
        File? models_table                       # TSV: model_id, model, trait_type, genetic_model, cat_var
        File? pheno_chunker_r                    # bin/pheno_chunker.R
        File? prepare_projects_py                # bin/prepare_projects.py
        Int pheno_chunk_size = 50                # Max phenotypes per regenie run
        Float missing_tolerance = 0.1            # Max proportion of missing samples

        ## ====================================================================
        ## MODE 3: PROJECTS ARRAY INPUTS
        ## ====================================================================
        ## Provide an array of ProjectConfig structs. Each element defines
        ## a separate project with its own phenotype/covariate files and settings.
        ## This mirrors the Nextflow projects_table TSV parameter.
        Array[ProjectConfig]? projects

        ## ====================================================================
        ## SHARED INPUTS (all modes)
        ## ====================================================================
        String genotypes_build                   # hg19 or hg38
        Array[String] chromosomes                # ["1","2",...,"22"]
        String project_date = ""
        Boolean phenotypes_delete_missings = false
        Int maxCatLevels = 10
        String regenie_gwas_test = "additive"    # Default genetic model (single mode)

        # Genotyped array data for Step 1
        File genotypes_array_bed
        File genotypes_array_bim
        File genotypes_array_fam

        # Imputed data for GWAS (optional)
        Array[GenotypeFiles]? genotypes_imputed
        String genotypes_imputed_format = "bgen"

        # Rare variant data (optional)
        Array[GenotypeFiles]? genotypes_rarevar
        String genotypes_rarevar_format = "bgen"

        # LD panel for clumping
        String ld_panel = "NO_LD_FILE"
        Array[File]? ld_panel_beds
        Array[File]? ld_panel_bims
        Array[File]? ld_panel_fams
        Array[String]? ld_panel_chroms

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
        String genes_group = "protein_coding"

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
        Int step2_gwas_chunk_size = 100000
        File? make_chunks_sql

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

        # Report settings
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
    ## SHARED: Compute GWAS variant chunks (once, shared across all projects)
    ## ========================================================================
    ## Note: rare variant chunks are computed per-project inside RunSingleAnalysis
    ## because MakeGenesChunks outputs File arrays which need File-type transfer.

    if (defined(genotypes_imputed) && step2_gwas_split && defined(make_chunks_sql)) {
        Array[GenotypeFiles] gwas_geno_for_chunks = select_first([genotypes_imputed])
        scatter (geno in gwas_geno_for_chunks) {
            call Tasks.MakeVariantsChunks {
                input:
                    filename = geno.file_prefix,
                    primary_file = geno.primary_file,
                    secondary_file = geno.secondary_file,
                    tertiary_file = geno.tertiary_file,
                    chromosome = geno.chromosome,
                    snplist_file = geno.secondary_file,
                    make_chunks_sql = select_first([make_chunks_sql]),
                    step2_gwas_chunk_size = step2_gwas_chunk_size,
                    chromosomes = chromosomes,
                    snplist_type = genotypes_imputed_format
            }
        }
        # Read chunk coordinates for each genotype file
        scatter (chunk_idx in range(length(gwas_geno_for_chunks))) {
            Array[String] computed_gwas_chunks_per_file = read_lines(MakeVariantsChunks.chunks_file[chunk_idx])
        }
    }

    ## ========================================================================
    ## MODE 1: SINGLE PROJECT
    ## ========================================================================
    ## Active when neither models_table nor projects is set.

    if (!defined(models_table) && !defined(projects)) {
        call RunAnalysisWF.RunSingleAnalysis as SingleRun {
            input:
                project_id = select_first([project]),
                phenotypes_file = select_first([phenotypes_filename]),
                phenotypes_columns = select_first([phenotypes_columns]),
                phenotypes_binary_trait = select_first([phenotypes_binary_trait]),
                regenie_gwas_test = regenie_gwas_test,
                covariates_file = covariates_filename,
                covariates_columns = covariates_columns,
                covariates_cat_columns = covariates_cat_columns,
                interaction_cov = interaction_cov,
                interaction_snp = interaction_snp,
                condition_list = condition_list,
                extract_snps_list = extract_snps_list,
                extract_genes_list = extract_genes_list,
                genotypes_array_bed = genotypes_array_bed,
                genotypes_array_bim = genotypes_array_bim,
                genotypes_array_fam = genotypes_array_fam,
                genotypes_imputed = genotypes_imputed,
                genotypes_imputed_format = genotypes_imputed_format,
                genotypes_rarevar = genotypes_rarevar,
                genotypes_rarevar_format = genotypes_rarevar_format,
                gwas_chunks_per_geno = computed_gwas_chunks_per_file,
                ld_panel = ld_panel,
                ld_panel_beds = ld_panel_beds,
                ld_panel_bims = ld_panel_bims,
                ld_panel_fams = ld_panel_fams,
                ld_panel_chroms = ld_panel_chroms,
                genotypes_build = genotypes_build,
                chromosomes = chromosomes,
                project_date = project_date,
                phenotypes_delete_missings = phenotypes_delete_missings,
                maxCatLevels = maxCatLevels,
                additional_geno_bgen_pgen_bed = additional_geno_bgen_pgen_bed,
                additional_geno_sample_psam_fam = additional_geno_sample_psam_fam,
                additional_geno_bgi_pvar_bim = additional_geno_bgi_pvar_bim,
                additional_geno_format = additional_geno_format,
                rarevar_set_list_file = rarevar_set_list_file,
                rarevar_anno_file = rarevar_anno_file,
                rarevar_mask_file = rarevar_mask_file,
                genes_bed = genes_bed,
                genes_ranges = genes_ranges,
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
                regenie_premade_predictions = regenie_premade_predictions,
                premade_prediction_files = premade_prediction_files,
                regenie_skip_predictions = regenie_skip_predictions,
                save_step1_predictions = save_step1_predictions,
                regenie_force_step1 = regenie_force_step1,
                regenie_ref_first_step1 = regenie_ref_first_step1,
                step1_use_loocv = step1_use_loocv,
                step1_niter = step1_niter,
                step1_n_chunks = step1_n_chunks,
                regenie_bsize_step2 = regenie_bsize_step2,
                regenie_ref_first_step2 = regenie_ref_first_step2,
                regenie_min_imputation_score = regenie_min_imputation_score,
                regenie_gwas_min_mac = regenie_gwas_min_mac,
                regenie_firth = regenie_firth,
                regenie_firth_approx = regenie_firth_approx,
                regenie_range = regenie_range,
                step2_gwas_split = step2_gwas_split,
                step2_rarevar_split = step2_rarevar_split,
                step2_rarevar_chunk_size = step2_rarevar_chunk_size,
                regenie_rarevar_min_mac = regenie_rarevar_min_mac,
                rarevars_aaf_bins = rarevars_aaf_bins,
                rarevars_vc_test = rarevars_vc_test,
                rarevars_joint_test = rarevars_joint_test,
                rarevars_vc_maxAAF = rarevars_vc_maxAAF,
                regenie_build_mask = regenie_build_mask,
                rarevars_write_mask_snplist = rarevars_write_mask_snplist,
                annotation_min_log10p = annotation_min_log10p,
                annotation_interval_kb = annotation_interval_kb,
                clumping = clumping,
                clump_p1 = clump_p1,
                clump_p2 = clump_p2,
                clump_kb = clump_kb,
                rarevar_min_log10p = rarevar_min_log10p,
                rarevar_stat_test = rarevar_stat_test,
                rarevar_stat_test_threshold = rarevar_stat_test_threshold,
                make_report = make_report,
                manhattan_annotations = manhattan_annotations,
                regional_plot_window_kb = regional_plot_window_kb,
                n_top_loci_plot = n_top_loci_plot,
                gwas_report_template = gwas_report_template,
                rarevar_report_template = rarevar_report_template,
                quarto_report_css = quarto_report_css,
                collapse_closegenes_py = collapse_closegenes_py,
                pipeline_version = pipeline_version
        }
    }

    ## ========================================================================
    ## MODE 2: MODELS TABLE
    ## ========================================================================
    ## Active when models_table is set.
    ## Runs pheno_chunker.R + prepare_projects.py to generate per-model
    ## phenotype/covariate files, then scatters over models.

    if (defined(models_table)) {
        call Tasks.SetupMultipleRuns {
            input:
                pheno_chunker_r = select_first([pheno_chunker_r]),
                prepare_projects_py = select_first([prepare_projects_py]),
                traits_table = select_first([phenotypes_filename]),
                models_table = select_first([models_table]),
                fam_file = genotypes_array_fam,
                pheno_chunk_size = pheno_chunk_size,
                missing_tolerance = missing_tolerance
        }

        scatter (models_proj_idx in range(length(SetupMultipleRuns.project_ids))) {
            # Conditionally set covariate file (File?) based on has_cov flag
            Boolean models_has_cov = (SetupMultipleRuns.has_cov_list[models_proj_idx] == "true")
            if (models_has_cov) {
                File models_real_cov_file = SetupMultipleRuns.cov_files[models_proj_idx]
            }

            call RunAnalysisWF.RunSingleAnalysis as ModelsRun {
                input:
                    project_id = SetupMultipleRuns.project_ids[models_proj_idx],
                    phenotypes_file = SetupMultipleRuns.pheno_files[models_proj_idx],
                    phenotypes_columns = SetupMultipleRuns.pheno_cols_list[models_proj_idx],
                    phenotypes_binary_trait = (SetupMultipleRuns.pheno_binary_list[models_proj_idx] == "True"),
                    regenie_gwas_test = SetupMultipleRuns.pheno_model_list[models_proj_idx],
                    covariates_file = models_real_cov_file,
                    covariates_columns = SetupMultipleRuns.cov_cols_list[models_proj_idx],
                    covariates_cat_columns = SetupMultipleRuns.cov_cat_cols_list[models_proj_idx],
                    interaction_cov = interaction_cov,
                    interaction_snp = interaction_snp,
                    condition_list = condition_list,
                    extract_snps_list = extract_snps_list,
                    extract_genes_list = extract_genes_list,
                    genotypes_array_bed = genotypes_array_bed,
                    genotypes_array_bim = genotypes_array_bim,
                    genotypes_array_fam = genotypes_array_fam,
                    genotypes_imputed = genotypes_imputed,
                    genotypes_imputed_format = genotypes_imputed_format,
                    genotypes_rarevar = genotypes_rarevar,
                    genotypes_rarevar_format = genotypes_rarevar_format,
                    gwas_chunks_per_geno = computed_gwas_chunks_per_file,
                    ld_panel = ld_panel,
                    ld_panel_beds = ld_panel_beds,
                    ld_panel_bims = ld_panel_bims,
                    ld_panel_fams = ld_panel_fams,
                    ld_panel_chroms = ld_panel_chroms,
                    genotypes_build = genotypes_build,
                    chromosomes = chromosomes,
                    project_date = project_date,
                    phenotypes_delete_missings = phenotypes_delete_missings,
                    maxCatLevels = maxCatLevels,
                    additional_geno_bgen_pgen_bed = additional_geno_bgen_pgen_bed,
                    additional_geno_sample_psam_fam = additional_geno_sample_psam_fam,
                    additional_geno_bgi_pvar_bim = additional_geno_bgi_pvar_bim,
                    additional_geno_format = additional_geno_format,
                    rarevar_set_list_file = rarevar_set_list_file,
                    rarevar_anno_file = rarevar_anno_file,
                    rarevar_mask_file = rarevar_mask_file,
                    genes_bed = genes_bed,
                    genes_ranges = genes_ranges,
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
                    regenie_premade_predictions = regenie_premade_predictions,
                    premade_prediction_files = premade_prediction_files,
                    regenie_skip_predictions = regenie_skip_predictions,
                    save_step1_predictions = save_step1_predictions,
                    regenie_force_step1 = regenie_force_step1,
                    regenie_ref_first_step1 = regenie_ref_first_step1,
                    step1_use_loocv = step1_use_loocv,
                    step1_niter = step1_niter,
                    step1_n_chunks = step1_n_chunks,
                    regenie_bsize_step2 = regenie_bsize_step2,
                    regenie_ref_first_step2 = regenie_ref_first_step2,
                    regenie_min_imputation_score = regenie_min_imputation_score,
                    regenie_gwas_min_mac = regenie_gwas_min_mac,
                    regenie_firth = regenie_firth,
                    regenie_firth_approx = regenie_firth_approx,
                    regenie_range = regenie_range,
                    step2_gwas_split = step2_gwas_split,
                    step2_rarevar_split = step2_rarevar_split,
                    step2_rarevar_chunk_size = step2_rarevar_chunk_size,
                    regenie_rarevar_min_mac = regenie_rarevar_min_mac,
                    rarevars_aaf_bins = rarevars_aaf_bins,
                    rarevars_vc_test = rarevars_vc_test,
                    rarevars_joint_test = rarevars_joint_test,
                    rarevars_vc_maxAAF = rarevars_vc_maxAAF,
                    regenie_build_mask = regenie_build_mask,
                    rarevars_write_mask_snplist = rarevars_write_mask_snplist,
                    annotation_min_log10p = annotation_min_log10p,
                    annotation_interval_kb = annotation_interval_kb,
                    clumping = clumping,
                    clump_p1 = clump_p1,
                    clump_p2 = clump_p2,
                    clump_kb = clump_kb,
                    rarevar_min_log10p = rarevar_min_log10p,
                    rarevar_stat_test = rarevar_stat_test,
                    rarevar_stat_test_threshold = rarevar_stat_test_threshold,
                    make_report = make_report,
                    manhattan_annotations = manhattan_annotations,
                    regional_plot_window_kb = regional_plot_window_kb,
                    n_top_loci_plot = n_top_loci_plot,
                    gwas_report_template = gwas_report_template,
                    rarevar_report_template = rarevar_report_template,
                    quarto_report_css = quarto_report_css,
                    collapse_closegenes_py = collapse_closegenes_py,
                    pipeline_version = pipeline_version
            }
        }
    }

    ## ========================================================================
    ## MODE 3: PROJECTS ARRAY
    ## ========================================================================
    ## Active when projects is set (and models_table is NOT set).
    ## Each ProjectConfig in the array defines a separate project.

    if (!defined(models_table) && defined(projects)) {
        Array[ProjectConfig] projects_list = select_first([projects])

        scatter (proj in projects_list) {
            call RunAnalysisWF.RunSingleAnalysis as ProjectsRun {
                input:
                    project_id = proj.project_id,
                    phenotypes_file = proj.phenotype_file,
                    phenotypes_columns = proj.phenotype_columns,
                    phenotypes_binary_trait = proj.phenotype_binary,
                    regenie_gwas_test = proj.genetic_model,
                    covariates_file = proj.covariate_file,
                    covariates_columns = proj.covariate_columns,
                    covariates_cat_columns = proj.covariate_cat_columns,
                    interaction_cov = proj.interaction_cov,
                    interaction_snp = proj.interaction_snp,
                    condition_list = proj.condition_list,
                    extract_snps_list = proj.extract_snps_list,
                    extract_genes_list = proj.extract_genes_list,
                    genotypes_array_bed = genotypes_array_bed,
                    genotypes_array_bim = genotypes_array_bim,
                    genotypes_array_fam = genotypes_array_fam,
                    genotypes_imputed = genotypes_imputed,
                    genotypes_imputed_format = genotypes_imputed_format,
                    genotypes_rarevar = genotypes_rarevar,
                    genotypes_rarevar_format = genotypes_rarevar_format,
                    gwas_chunks_per_geno = computed_gwas_chunks_per_file,
                    ld_panel = ld_panel,
                    ld_panel_beds = ld_panel_beds,
                    ld_panel_bims = ld_panel_bims,
                    ld_panel_fams = ld_panel_fams,
                    ld_panel_chroms = ld_panel_chroms,
                    genotypes_build = genotypes_build,
                    chromosomes = chromosomes,
                    project_date = project_date,
                    phenotypes_delete_missings = phenotypes_delete_missings,
                    maxCatLevels = maxCatLevels,
                    additional_geno_bgen_pgen_bed = additional_geno_bgen_pgen_bed,
                    additional_geno_sample_psam_fam = additional_geno_sample_psam_fam,
                    additional_geno_bgi_pvar_bim = additional_geno_bgi_pvar_bim,
                    additional_geno_format = additional_geno_format,
                    rarevar_set_list_file = rarevar_set_list_file,
                    rarevar_anno_file = rarevar_anno_file,
                    rarevar_mask_file = rarevar_mask_file,
                    genes_bed = genes_bed,
                    genes_ranges = genes_ranges,
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
                    regenie_premade_predictions = regenie_premade_predictions,
                    premade_prediction_files = premade_prediction_files,
                    regenie_skip_predictions = regenie_skip_predictions,
                    save_step1_predictions = save_step1_predictions,
                    regenie_force_step1 = regenie_force_step1,
                    regenie_ref_first_step1 = regenie_ref_first_step1,
                    step1_use_loocv = step1_use_loocv,
                    step1_niter = step1_niter,
                    step1_n_chunks = step1_n_chunks,
                    regenie_bsize_step2 = regenie_bsize_step2,
                    regenie_ref_first_step2 = regenie_ref_first_step2,
                    regenie_min_imputation_score = regenie_min_imputation_score,
                    regenie_gwas_min_mac = regenie_gwas_min_mac,
                    regenie_firth = regenie_firth,
                    regenie_firth_approx = regenie_firth_approx,
                    regenie_range = regenie_range,
                    step2_gwas_split = step2_gwas_split,
                    step2_rarevar_split = step2_rarevar_split,
                    step2_rarevar_chunk_size = step2_rarevar_chunk_size,
                    regenie_rarevar_min_mac = regenie_rarevar_min_mac,
                    rarevars_aaf_bins = rarevars_aaf_bins,
                    rarevars_vc_test = rarevars_vc_test,
                    rarevars_joint_test = rarevars_joint_test,
                    rarevars_vc_maxAAF = rarevars_vc_maxAAF,
                    regenie_build_mask = regenie_build_mask,
                    rarevars_write_mask_snplist = rarevars_write_mask_snplist,
                    annotation_min_log10p = annotation_min_log10p,
                    annotation_interval_kb = annotation_interval_kb,
                    clumping = clumping,
                    clump_p1 = clump_p1,
                    clump_p2 = clump_p2,
                    clump_kb = clump_kb,
                    rarevar_min_log10p = rarevar_min_log10p,
                    rarevar_stat_test = rarevar_stat_test,
                    rarevar_stat_test_threshold = rarevar_stat_test_threshold,
                    make_report = make_report,
                    manhattan_annotations = manhattan_annotations,
                    regional_plot_window_kb = regional_plot_window_kb,
                    n_top_loci_plot = n_top_loci_plot,
                    gwas_report_template = gwas_report_template,
                    rarevar_report_template = rarevar_report_template,
                    quarto_report_css = quarto_report_css,
                    collapse_closegenes_py = collapse_closegenes_py,
                    pipeline_version = pipeline_version
            }
        }
    }

    ## ========================================================================
    ## OUTPUTS
    ## ========================================================================
    ## Outputs are organized by mode. Only the active mode produces results.

    output {
        ## ---- Single project mode outputs ----
        Array[File]? single_step1_predictions = SingleRun.step1_predictions
        File? single_step1_log = SingleRun.step1_log
        File? single_validated_phenotypes = SingleRun.validated_phenotypes
        File? single_phenotype_validation_log = SingleRun.phenotype_validation_log
        File? single_covariate_validation_log = SingleRun.covariate_validation_log
        Array[File]? single_gwas_merged_results = SingleRun.gwas_merged_results
        Array[File]? single_gwas_merged_results_tbi = SingleRun.gwas_merged_results_tbi
        File? single_gwas_step2_log = SingleRun.gwas_step2_log
        Array[File]? single_gwas_filtered_results = SingleRun.gwas_filtered_results
        Array[File]? single_gwas_annotated_results = SingleRun.gwas_annotated_results
        Array[File?]? single_gwas_toploci = SingleRun.gwas_toploci
        Array[File?]? single_gwas_annotloci = SingleRun.gwas_annotloci
        Array[File]? single_gwas_reports = SingleRun.gwas_reports
        Array[File]? single_rarevar_merged_results = SingleRun.rarevar_merged_results
        Array[File]? single_rarevar_merged_results_tbi = SingleRun.rarevar_merged_results_tbi
        File? single_rarevar_step2_log = SingleRun.rarevar_step2_log
        Array[File]? single_rarevar_filtered_results = SingleRun.rarevar_filtered_results
        Array[File]? single_rarevar_processed_results = SingleRun.rarevar_processed_results
        Array[File]? single_rarevar_reports = SingleRun.rarevar_reports

        ## ---- Models table mode outputs (one array element per model) ----
        File? models_analysis_config = SetupMultipleRuns.analysis_config
        File? models_master_table = SetupMultipleRuns.master_table
        Array[Array[File]]? models_step1_predictions = ModelsRun.step1_predictions
        Array[File?]? models_step1_logs = ModelsRun.step1_log
        Array[File]? models_validated_phenotypes = ModelsRun.validated_phenotypes
        Array[File]? models_phenotype_validation_logs = ModelsRun.phenotype_validation_log
        Array[File?]? models_covariate_validation_logs = ModelsRun.covariate_validation_log
        Array[Array[File]?]? models_gwas_merged_results = ModelsRun.gwas_merged_results
        Array[Array[File]?]? models_gwas_merged_results_tbi = ModelsRun.gwas_merged_results_tbi
        Array[File?]? models_gwas_step2_logs = ModelsRun.gwas_step2_log
        Array[Array[File]?]? models_gwas_filtered_results = ModelsRun.gwas_filtered_results
        Array[Array[File]?]? models_gwas_annotated_results = ModelsRun.gwas_annotated_results
        Array[Array[File]?]? models_gwas_reports = ModelsRun.gwas_reports
        Array[Array[File]?]? models_rarevar_merged_results = ModelsRun.rarevar_merged_results
        Array[File?]? models_rarevar_step2_logs = ModelsRun.rarevar_step2_log
        Array[Array[File]?]? models_rarevar_filtered_results = ModelsRun.rarevar_filtered_results
        Array[Array[File]?]? models_rarevar_reports = ModelsRun.rarevar_reports

        ## ---- Projects array mode outputs (one array element per project) ----
        Array[Array[File]]? projects_step1_predictions = ProjectsRun.step1_predictions
        Array[File?]? projects_step1_logs = ProjectsRun.step1_log
        Array[File]? projects_validated_phenotypes = ProjectsRun.validated_phenotypes
        Array[File]? projects_phenotype_validation_logs = ProjectsRun.phenotype_validation_log
        Array[File?]? projects_covariate_validation_logs = ProjectsRun.covariate_validation_log
        Array[Array[File]?]? projects_gwas_merged_results = ProjectsRun.gwas_merged_results
        Array[Array[File]?]? projects_gwas_merged_results_tbi = ProjectsRun.gwas_merged_results_tbi
        Array[File?]? projects_gwas_step2_logs = ProjectsRun.gwas_step2_log
        Array[Array[File]?]? projects_gwas_filtered_results = ProjectsRun.gwas_filtered_results
        Array[Array[File]?]? projects_gwas_annotated_results = ProjectsRun.gwas_annotated_results
        Array[Array[File]?]? projects_gwas_reports = ProjectsRun.gwas_reports
        Array[Array[File]?]? projects_rarevar_merged_results = ProjectsRun.rarevar_merged_results
        Array[File?]? projects_rarevar_step2_logs = ProjectsRun.rarevar_step2_log
        Array[Array[File]?]? projects_rarevar_filtered_results = ProjectsRun.rarevar_filtered_results
        Array[Array[File]?]? projects_rarevar_reports = ProjectsRun.rarevar_reports
    }
}
