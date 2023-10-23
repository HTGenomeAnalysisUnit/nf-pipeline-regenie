# Full list of parameters

## Global settings of resources

Global settings for max resources available per process in your system

- max_memory 
- max_cpus   
- max_time   

## General settings per run

- project: project name. **Max 50 chars. No spaces allowed. Use only alphanumeric, underscores, or dashes**
- project_date: usually configured automatically to execution date, but you can override the date manually
- outdir: a specific outdir for the project, if null the project id is used as output dir, otherwise `outdir/project_id/`
- chromosomes: chromosomes to include in the analysis. Accepts comma separated lists and intervals like `1-5,8-12,22`
- master_log_dir: set this to redirect all logs to a specific folder

## Multi-model specific settings
    
	models_table                          = null //when a table is provided, this activates master table mode
	pheno_chunk_size                      = 50 //maximum number of phenotypes for a single GWAS run
	missing_tolerance                     = 0.1 //max fraction of missing values per phenotype when assembling execution groups

## Input files

### Genetic data

    genotypes_array                       = null //genotypes data for step 1. Can be bgen, pgen, bed or vcf.gz
    genotypes_imputed                     = null //genotype data for gwas analysis. Can be bgen, pgen, bed or vcf.gz
    genotypes_rarevar                     = null //genotype data for rare variants analysis. Can be bgen, pgen, plink bed or vcf.gz
    genotypes_build                       = null //genome build
    genotypes_imputed_format              = null //input data format for gwas. Can be bgen, pgen, bed or vcf
    genotypes_rarevar_format              = null //input data format for rare variants. Can be bgen, pgen, bed or vcf
    imputed_sample_file                   = 'NO_SAMPLE_FILE' //Provide a specific sample file to use with gwas bgen input
    rarevar_sample_file                   = 'NO_SAMPLE_FILE' //Provide a specific sample file to use with rare variants bgen input
    ld_panel                              = 'NO_LD_FILE' //A pattern pointing to a subset of genotypes to be used for LD computation. Optional, but highly reccomended for large datasets

### Phenotypes

    phenotypes_filename                   = null //table of phenos - required
    phenotypes_columns                    = null //comma separated list of col names in pheno file to be analyzed
    phenotypes_binary_trait               = null //true for binary
    
### Covariates

    covariates_filename                   = 'NO_COV_FILE' //files containing covariates, use NO_COV_FILE when absent
    covariates_columns                    = '' //comma-separated list of covariates column names 
    covariates_cat_columns                = '' //comma-separated list of categorical covariate column names
    maxCatLevels                          = 10 //maximum number of allowed levels for categorical covars

### Rare variants accessory files
These are mandatory when running rare var analysis

    rarevar_set_list_file                 = null //set list file as defined in regenie docs
    rarevar_anno_file                     = null //variant annotation file as defined in regenie docs
    rarevar_mask_file                     = null //mask definition file as defined in regenie docs

### Other optional inputs
    genes_bed                             = false //an optional .bed file specifying genes intervals used for results annotation
    genes_ranges                          = false //an optional .interval file specifying genes intervals used for loci annotation

## STEP1 PRE-PROCESSING SETTINGS

### SNP pruning

    prune_enabled                         = false
    prune_maf                             = 0.01
    prune_window_kbsize                   = 1000
    prune_step_size                       = 100
    prune_r2_threshold                    = 0.9

### Filtering

    qc_maf                                = '0.01'
    qc_mac                                = '100'
    qc_geno                               = '0.05'
    qc_hwe                                = '1e-15'
    qc_mind                               = '0.1'

## REGENIE STEP1 SETTINGS

    regenie_bsize_step1                   = 1000
    regenie_premade_predictions           = false //or pattern to regenie step1 files
    save_step1_predictions                = true
    regenie_force_step1                   = false
    regenie_ref_first_step1               = false
    step1_use_loocv                       = false
    step1_niter                           = 30
    step1_n_chunks                        = 100 // N chunks when performing step1 L0 regression
    
## REGENIE STEP2 SETTINGS

### General settings

    phenotypes_delete_missings            = false //remove samples with missing data at any of the phenotypes
    regenie_bsize_step2                   = 400
    regenie_ref_first_step2               = true
    regenie_skip_predictions              = false //skip reading the step1 predictions (corresponds to simple linear/logistic regression)
    regenie_range                         = '' // when splitting is not active you can use this to specify a genomic range for step2 analysis
    regenie_extract_snps                  = '' // when splitting is not active you can specify a file containing a list of variant IDs to restrict step2 analysis
    regenie_extract_genes                 = '' // when splitting is not active you can use this to specify a file containing a list of genes to restrict step2 analysis
    interaction_cov                       = null // run GxE test in GWAS specifying the interacting covariate from covariate table
    interaction_snp	                      = null // run GxG test in GWAS specifying the interacting variant ID
    condition_list                        = null // run conditional analysis in GWAS specifying a files with variant IDs to condition on
    save_chunks_file                      = true
    save_step2_logs                       = true
    save_bgen_index                       = true
    save_bgen_sample                      = true    
    save_pgen                             = true

### VCF to PGEN conversion settings

    gwas_read_dosage_from                 = 'DS' //DS (usually for VCF from imputateion) or GP (usually VCF from sequencing)
    rarevar_read_dosage_from              = 'GP' //DS (usually for VCF from imputateion) or GP (usually VCF from sequencing)
    import_dosage_certainty               = 0.7 //when using GP, the certainty threshold to import dosages. If none of the probabilities is above this, the genotype is set to missing
    vcf_fixed_fid                         = null //when null vcf is converted to pgen using --double-id, otherwise fid is fixed to this value

### GWAS analysis settings

    step2_gwas_split                      = true //when true activate split of step2 by variant chunks
    step2_gwas_chunk_size                 = 100000 //n variants per chunk when running gwas in split mode
    regenie_gwas_test                     = 'additive' //additive, dominant, recessive
    regenie_min_imputation_score          = '0.00'
    regenie_gwas_min_mac                  = '50' // min MAC for variants to be included in step2 for gwas 
    regenie_firth                         = true
    regenie_firth_approx                  = true

### Rare variant analysis settings

    step2_rarevar_split                   = true //when true activate split of step2 by gene chunks
    step2_rarevar_chunk_size              = 200 //n genes per chunk when running rare variant test in split mode
    regenie_rarevar_min_mac               = '1' // min MAC for variants to be included in step2 for rare vars
    rarevars_aaf_bins                     = '0.01,0.05' //comma-separated list of AAF upper bounds to use when building masks for burden test
    rarevars_vc_test                      = 'skat,skato,acatv,acato' //comma-separated list of SKAT/ACAT-type tests to run
    rarevars_vc_maxAAF                    = '0.05' //AAF upper bound to use for SKAT/ACAT-type tests
    regenie_build_mask                    = 'max' //build mask for rare variant test. Can be max, sum, or a comphet
    rarevars_write_mask_snplist           = false //when true write list of variants that went into each mask to file

## GWAS RESULTS ANNOTATION AND CLUMPING

    genes_group                           = 'protein_coding' //genes group to use for annotation. Can be all or protein_coding
    annotation_min_log10p                 = 7.3 //results with -log10(p) above this will be reported as top hits with annotated genes
    annotation_interval_kb                = 25
    clumping                              = true
    clump_p1                              = 5e-8
    clump_p2                              = 1e-4
    clump_kb                              = 250
    
## RARE VARIANTS RESULTS ANNOTATION

    rarevar_min_log10p                    = 5 //results with -log10(p) above this will be reported as top hits
    rarevar_stat_test                     = "BONF_bygroup" //Stat value to filter on. Possible values: "FDR_bygroup", "FDR_alltests", "BONF_bygroup", "BONF_alltests"
    rarevar_stat_test_threshold           = 1.3 //results with -log10(stat_value) above this will be annotated in a dedicated manhattan plot

## REPORT SETTINGS

    make_report                           = true //it can be useful to disable when analyzing many phenotypes
    manhattan_annotations                 = 'genes' //how to annotate peaks in manhattan plot. Either 'genes' or 'snpid'
    regional_plot_window_kb               = 300 //window size in kb for regional plots. This value is added on each side of the locus when plotting
    n_top_loci_plot                       = 5 //number of top loci to plot in regional plots
