# Phenotype and covariates

To perform the analysis you should provide one tab-separated file containing the phenotypes(s) of interest and optionally and additional tab-separated files containing the covariates to be included in the analysis. In the end, only samples found in both phenotype and covariates tables will be processed.

## Phenotypes

- The phenotype input table is defined with the `phenotypes_filename` parameter and must be a tab-separated file with header defining the phenotype(s). The first 2 columns must be FID, IID.
- The exact list of column names to process are then defined with the `phenotypes_columns` parameters. 
- Set `phenotypes_binary_trait` to true or false according to the type of phenotype(s). Note that quantitative and binary traits can not be mixed in single mode run.

### Notes for multi-models run

- In multi-models mode, all phenotypes and covariates values must be provided in a single tab-separated file with header and configured as `phenotypes_filename`. The exact data types and column names for phenotypes and covariates are then configured in the models table []().
- When using multi-model mode, the phenotype file can have a single `IID` column as first column. In this case, the `FID` will be deduced from the fam file in the `genotypes_array` dataset.

## Covariates

- The covariates table is defined with the `covariates_filename` parameter and must be a tab-separated file with header defining the covariate(s). You can set this to 'NO_COVAR_FILE' when no covariates are present in your model. 
- The exact list of column names to process are then defined with the `covariates_columns`. Here it is fine to mix binary and quantitative covars. 
- Binary covars must contain only 0/1 values (otherwise they will be treated as numeric values).
- Columns containing categorical values must be specified with the `covariates_cat_columns` parameter. These will be automatically converted to dummy variables. The maximum number of categories allowed is 10 by default, but can be adjusted with the `maxCatLevels` parameter.

Note that **no missing values** are allowed in this table.

## Processing of missing values

- Samples listed in pheno file that are not in bgen/bed/pgen file are ignored. Genotyped samples that are not in the pheno file are removed from the analysis.
- With quantitative traits, missing values are mean-imputed in Step 1 and they are dropped when testing each phenotype in Step 2
- With binary traits, missing values are mean-imputed in Step 1 when fitting the level 0 linear ridge regression and they are dropped when fitting the level 1 logistic ridge regression for each trait. In Step 2, missing values are dropped when testing each trait.