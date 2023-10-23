# Projects table

When executing a multi-projects run you need to configure the genetic input data in the main configuration file, the other project-specific parameters will be taken from the projects table.

The parameter `projects_table` is used to activate a multi-projects run by configuring a projects table. The projects table is a tab-separated file with header and the following columns (names are mandatory):

- `project_id`: a unique identifier for the project
- `pheno_file`: path to the phenotype file
- `pheno_cols`: comma-separated list of phenotype columns
- `pheno_binary`: comma-separated list of binary phenotype columns
- `pheno_model`: comma-separated list of phenotype models
- `cov_file`: path to the covariates file, use NO_COV_FILE if no covariates are used
- `cov_cols`: comma-separated list of covariates columns (excluding categorical covariates), can be omitted if cov_file is NO_COV_FILE
- `cov_cat_cols` (optional): comma-separated list of categorical covariates columns
- `interaction_snp` (optional): a single variant ID to perform interaction analysis (GxG)
- `interaction_cov` (optional): a single covariate name to perform interaction analysis (GxE)
- `condition_list` (optional): a text file containing a list of variants IDs to condition on

An example of a models table is shown below:

```bash
model_id        model   trait_type      genetic_model cat_var
M1      QP1 ~ Q1+Q2+Q5  quant   additive  NA
M2      QP2 ~ Q1+Q2+Q5  quant   additive  Q1,Q2
M3      QP3 ~ Q1+Q2+Q5  quant   additive  NA
M4      QP4 ~ 1  quant   additive   NA
```
