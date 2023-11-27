# Main parameters

A template for preparing a configuration file for a new project is provided in `templates/single_project.conf`. You can copy this file and edit it to adjust the parameters for your analysis.

Here is a list of the main parameters you need to adjust for a new analysis:

- Set `project` to a custom string defining your project ID. Note that max 50 chars are allowed and no special characters. This will be used in reports and a folder will be created with the same name to store all results.

- Set `chromosomes` to represent the list of chromosome to be included in the analysis. You can use a comma-separated list of chromosome numbers, a range like `1-22` or a mix like `1,4,11-18`.

- You can eventually restrict the analysis to a specific genomic range, or a specific list of SNP IDs or gene IDs using the `regenie_range`, `regenie_extract_snps` and `regenie_extract_genes` options. This requires `step2_gwas_split` and `step2_rarevar_split` to be set to `false`. Keep in mind that in this case you must ensure that the provided SNPs/genes/region are present in the genotype dataset you provided as input and in the specified chromosomes. Otherwise, the pipeline will fail due to SNPs remaining for step2 analysis.

- Set `genotypes_build` to the build of your genotype data, either hg19 or hg38.

- Set `step2_gwas_chunk_size` and `step2_rarevar_chunk_size` according to the size of your dataset. These control how the dataset is split for step2 analysis. Teh default values usually work fine, but you can increase/decrease them if you have very large or small datasets. Keep in mind that a small chunk size may result in a very large amount of parallel jobs. By default the pipeline job submission rate is limited to 200 concurrent jobs. The total number of jobs will be *N_SNPs/gwas_chunk_size* for GWAS and *N_genes/rare_chunk_size* for rare variant analysis.

- For GWAS analysis, set `regenie_gwas_test` to the type of model you want to test, either 'additive', 'dominant' or 'recessive'.
  
- For rare variant test, set `rarevars_vc_test` as a comma-separated list of rare variants tests to perform among those accepted by regenie: skat, skato, skato-acat, acatv, acato, acato-full.

- Set `regenie_gwas_min_mac` and `regenie_rarevar_min_mac` to control the min allowed MAC for variants tested in either GWAS or rare variant analysis. Variants with MAC below this threshold will be excluded from the analysis.

- Set `annotation_min_log10p` to the min value ( -log10(pval) ) for top hit SNPs from GWAS results. These SNPs are also annotated in the manhattan plot in the HTML report.

- Set `rarevar_min_log10p` to the min value ( -log10(pval) ) for top hit genes from rare variants analysis. These genes are also annotated in the general manhattan plot in the HTML report.

- Set `clump_p1` to the maximum pvalue allowed for index SNPs during plink clumping to define top loci

- If you are analyzing many phenotypes on a large dataset, we suggest to set `make_report` to false or at least disable the generation of the locus zoom plots by setting `n_top_loci_plot` to zero. Generating the HTML graphical report and especially the regional plots, can add a considerable amount of time when many phenotypes are tested together on a large dataset.

- If you have categorical covariates, the maximum number of allowed categories is set to 10 bu default. You can adjsut this using the `maxCatLevels` parameter.

## Input files

Set the input files parameters as described in the input files sections.

- The [full genetic dataset](input-full-data.md) to perform GWAS analysis
- The [reduced SNPs dataset](input-indep-snps.md) to perform REGENIE step1 regression
- The [rare variants dataset and annotations](input-rarevars-data.md) to perform rare variant tests
- An optional set of files representing your [LD panel](input-ld-panel.md) to perform loci clumping

## Multi-models run

In case you want to configure a you also need to set the following parameters:

- Set `models_table` to a tab-separated file defining the models to test
  
- Set `missing_tolerance` to the maximum allowed fraction of missing phenotype values when collecting uniform group of phenotypes for a run.

## Multi-projects run

In case you want to configure a multi-projects run you also need to set the following parameters:

- Set `projects_table` to a tab-separated file defining the projects to test

## Conditional and interaction analysis

In case you want to perform conditional or interaction analysis, you can configure the following parameters:

- `interaction_cov`: run GxE test specifying the interacting covariate from covariate table
- `interaction_snp`: run GxG test specifying the interacting variant ID
- `condition_list`: run conditional analysis specifying a file with variant IDs to condition on

Note that to perform conditional/interaction analysis, an additional genotype dataset must be provided containing the SNP(s) used for conditioning/interaction. This can be configured using the following parameters:

- `additional_geno_file`: prefix of the genotype dataset containing vars in condition_list or interaction var. This is mandatory for conditional or interaction analysis
- `additional_geno_format`: can be bgen, pgen or bed.

See the [conditional analysis section](input-conditional-analysis.md) for more details.
