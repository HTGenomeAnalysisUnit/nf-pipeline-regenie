# Single model mode

The single model execution mode represent a standard regenie assoiation analysis. Given a set of phenote variables all of the same type and a set of covariates you can perform both association test at single variant level (GWAS) or rare variants aggregated tests.

Ideally, none of the tested phenotype should have a large number of missing values (10-15% is  areasonable threshold) since missing values are automatically imputedby regenie in some steps and this may introduce biases.

You can test multiple phenotypes together in the same execution as long as they are of the same type, namely all binary or quantitativa variables, and they are all tested with the same set of covariates. Regenie is able to manage even alarge phenotypes batch, and more than 100 quantitative phenotypes can be tested together. Computation is more time demanding for bianry traits, so we sggested to reduce the phenotype batches in this scenario.

Given an input dataset containing imputed genotypes and/or sequencing data, the pipeline automatically split the analysis in chunks to maximize performances. Results are then concatenated by phenotype, filtered and annotated with nearby genes

All configuration files are saved in the results folder so that the analysis is easy to reproduce.
