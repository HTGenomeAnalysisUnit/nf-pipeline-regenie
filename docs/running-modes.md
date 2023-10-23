# Execution modes

## Single project

The single model execution mode represent a standard regenie assoiation analysis. Given a set of phenotypes variables all of the same type (quantitative or binary) and a set of covariates you can perform both association test at single variant level (GWAS) or rare variants aggregated tests.

Ideally, none of the tested phenotype should have a large number of missing values (10-15% is a reasonable threshold) since missing values are automatically imputed by regenie in some steps and this may introduce biases. Also keep in mind that samples with missing values in any of the covariates will be excluded from the analysis.

You can test multiple phenotypes together in the same execution as long as they are of the same type, namely all binary or quantitativa variables, and they are all tested with the same set of covariates. Regenie is able to manage large phenotypes batches in a single execution. When working with biobank scale datasets (>200k individuals and >20M variants) our suggestion is to limit each analyses to 100 quantitative phenotypes or 50 binary phenotypes, but you can run larger batches if you have enough computational resources or a smaller datasets.

In single project mode, all parameters for the analysis are taken from the configuration file. See the [main parameters section](main-parameters.md) for more details on how to prepare the configuration file.

## Multi-models

When using multi-models execution mode one can test any number of arbitrary models (combinations of pheotypes and covariates) on a given cohort. In this case, user needs to provide teh following parameters:

- `models_table`: a models tables describing the models of interest
- `phenotypes_filename`: a tab-separated file containing all the traits (pheotypes and covarites)
- `missing_tolerance`: a threshold for the maximum fraction of missing values in phenotype variables
- `pheno_chunk_size`: a threshold for the maximum size of a phenotype batch

Based on this information we first perform a preparing step that configure uniform analysis runs by grouping together compatible models. The process works as follows:

1. Models with the same phenotyped data type, the same genetic model and the same set of covariates are grouped together
2. For each model a matrix is generated containing all the phenotype variables and covariates used in the model
3. Such model matrixes are then grouped to create multi-phenotype analysis in which none of the tested phenotype has more then the allowed fraction of missing values and no missing values are present in any covariates
4. Analysis chunks are then prepared based on the configured phenotype batch size
5. The resulting analysis chunks are automatically submitted to the pipeline for association analysis

A main result folder is created under `outdir` (or using the default nf-regenie-pipeline_results) where models configuration is saved. Then results for each analysis chunk are saved in a separate sub-folder containing the configuration and the exact tables of phenotypes and covariates used in the execution.

In this way multiple hypotheses can be easily tested on large cohort with minimum amount of data pre-processing.

Note that in the output all phenotypes are renamed according to the model id defined in the input models table. In this way, each phenotype code in the output represents a unique model (combination of a phenotype and covariates). So for example, suppose you configured a multi-model run to test the phenotype `height` with 2 different sets of covariates:

| model_id | model |
| -------- | ----- |
| M1 | height ~ cov1+cov2 |
| M2 | height ~ cov3+cov4 |

The resulting output traits will be named `M1` and `M2` and the output files containing summary statistics will be named `M1.regenie.gz`, `M2.regenie.gz`.

See the [multi models configuration section](input-models-table.md) for more details on how to prepare the input files for a multi-models run.

## Multi-projects

Multi-projects mode can be used to conveniently perform multiple project runs on the same genetic data in a single execution. In this case, the user needs to provide the `projects_table` parameter to specify a tab-separated file describing the projects to run. See the [projects table section](input-projects-table.md) for more details on how to prepare this file.

Values from this table will be used to configure multiple analyses and will over-ride any value provided in the config files.

A main result folder is created under `outdir` (or using the default nf-regenie-pipeline_results) where models configuration is saved. Then results for each analysis chunk are saved in a separate sub-folder containing the configuration and the exact tables of phenotypes and covariates used in the execution.
