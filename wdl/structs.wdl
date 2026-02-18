version 1.0

## ============================================================================
## Structs for the REGENIE GWAS Pipeline (WDL conversion)
## ============================================================================

struct PhenoMeta {
    String cols          # comma-separated phenotype column names
    Boolean binary       # true for binary traits
    String model         # additive, dominant, recessive
}

struct CovarMeta {
    String cols          # comma-separated covariate column names
    String cat_cols      # comma-separated categorical covariate column names
    String? gxe          # GxE interaction covariate
    String? gxg          # GxG interaction variant ID
}

struct AccessoryFiles {
    File? condition_list
    File? additional_bgen_pgen_bed
    File? additional_sample_psam_fam
    File? additional_bgi_pvar_bim
    File? extract_snps_list
    File? extract_genes_list
}

struct ProjectData {
    String project_id
    File phenotype_file
    PhenoMeta pheno_meta
    File covariate_file
    CovarMeta covar_meta
    AccessoryFiles accessory_files
}

struct PlinkDataset {
    File bed
    File bim
    File fam
}

struct BgenDataset {
    File bgen
    File bgi
    File sample
}

struct GenotypeFiles {
    String file_prefix
    File primary_file      # bed, bgen, or pgen
    File secondary_file    # bim, bgi, or pvar
    File tertiary_file     # fam, sample, or psam
    String chromosome      # chromosome or "ONE_FILE"
}

struct RuntimeAttributes {
    Int cpu
    String memory
    Int disk_size_gb
    Int max_retries
    String docker
}

## ProjectConfig: Defines a single project for multi-project / multi-model runs.
## Used when providing an Array[ProjectConfig] in projects_table mode,
## or constructed internally by models_table mode.
struct ProjectConfig {
    String project_id
    File phenotype_file
    String phenotype_columns        # comma-separated phenotype column names
    Boolean phenotype_binary        # true for binary, false for quantitative
    String genetic_model            # additive, dominant, or recessive
    File? covariate_file            # omit or null if no covariates
    String covariate_columns        # comma-separated covariate columns (empty string = none)
    String covariate_cat_columns    # comma-separated categorical covariates (empty string = none)
    String? interaction_cov         # GxE interaction covariate name
    String? interaction_snp         # GxG interaction variant ID
    File? condition_list            # file with variant IDs to condition on
    File? extract_snps_list         # file with variant IDs to restrict GWAS to
    File? extract_genes_list        # file with gene IDs to restrict rare variant analysis to
}
