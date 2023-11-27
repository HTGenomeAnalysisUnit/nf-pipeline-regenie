# Conditional and interaction analysis

In case you want to perform conditional or interaction analysis, you can configure the following parameters:

- `interaction_cov`: run GxE test specifying the interacting covariate from covariate table.
  For GxE tests where the interacting variable is categorical, you can specify the baseline level using `'VARNAME[BASE_LEVEL]'` (e.g. `'BMI[<25]'`). Otherwise, the first value found in the covariate file will be used as the baseline level.

- `interaction_snp`: run GxG test specifying the interacting variant ID
  For GxG tests, the default coding for the interacting variant is additive. If you would like to use dominant/recessive/categorical coding, use `'SNP_NAME[dom/rec/cat]'` (for example with dominant coding, `'SNPNAME[dom]'` will allow for separate effects between carriers vs non-carriers of the interacting variant). The allowed values in the brackets are add/dom/rec/cat.

- `condition_list`: run conditional analysis specifying a file with variant IDs to condition on

Note that to perform conditional/interaction analysis, an additional genotype dataset must be provided containing the SNP(s) used for conditioning/interaction. This can be configured using the following parameters:

- `additional_geno_format`: can be bgen, pgen or bed.
- `additional_geno_file`: prefix of the genotype dataset containing vars in condition_list or interaction var. This is mandatory for conditional or interaction analysis. Depending on the value of `additional_geno_format`, the pipeline expects specific files to be present
  - bgen: `additional_geno_file.bgen`, `additional_geno_file.bgen.bgi`, `additional_geno_file.sample`
  - bed: `additional_geno_file.bed`, `additional_geno_file.bim`, `additional_geno_file.fam`
  - pgen: `additional_geno_file.pgen`, `additional_geno_file.psam`, `additional_geno_file.pvar`
