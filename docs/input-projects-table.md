# Projects table

When executing a multi-projects run you need to configure the genetic input data in the main configuration file, the other project-specific parameters will be taken from the projects table.

The parameter `projects_table` is used to activate a multi-projects run by configuring a projects table. The projects table is a tab-separated file with header and the following columns (names are mandatory). You can use NA in the optional columns when the value is not relevant.

- `project_id`: a **unique** identifier for the project
- `pheno_file`: path to the phenotype file
- `pheno_cols`: comma-separated list of phenotype columns
- `pheno_binary`: True if the phenotype is binary, False otherwise
- `pheno_model`: association model for this phenotype, can be additive, dominant or recessive
- `cov_file`: path to the covariates file, use NO_COV_FILE if no covariates are used
- `cov_cols`: comma-separated list of covariates columns (excluding categorical covariates), can be omitted if cov_file is NO_COV_FILE
- `cov_cat_cols` (optional): comma-separated list of categorical covariates columns
- `interaction_snp` (optional): a single variant ID to perform interaction analysis (GxG)
- `interaction_cov` (optional): a single covariate name to perform interaction analysis (GxE)
- `condition_list` (optional): a text file containing a list of variants IDs to condition on

## Examples

A minimal projects table can be as follows:

| project_id | pheno_file | pheno_cols | pheno_binary | pheno_model | cov_file | cov_cols |
|------------|------------|------------|--------------|-------------|----------|----------|
| project1   | phenos.txt | qpheno1,qpheno2     | False   | additive      | covars.txt | cov1,cov2     |
| project2   | phenos.txt | bpheno1,bpheno2     | True   | additive      | NO_COV_FILE | NA    |

An example project table with all columns can be as follows:

| project_id | pheno_file | pheno_cols | pheno_binary | pheno_model | cov_file | cov_cols | cov_cat_cols | interaction_snp | interaction_cov | condition_list |
|------------|------------|------------|--------------|-------------|----------|----------|--------------|-----------------|-----------------|----------------|
| project1   | phenos.txt | qpheno1,qpheno2     | False   | additive      | covars.txt | cov1,cov2     | cat_covar1      | NA              | NA              | NA             |
| project2   | phenos.txt | qpheno1,qpheno2     | False   | additive      | covars.txt | cov1,cov2     | cat_covar1      | rs12345             | NA              | conditional_snps.txt             |
| project3   | phenos.txt | bpheno3,bpheno4     | True   | additive    | NO_COV_FILE | NA     | NA     | NA              | NA              | NA             |
