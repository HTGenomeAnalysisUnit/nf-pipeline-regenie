# Multi-models exeution

When using multi-models eecution mode one can test any number of arbitrary models (cmbinations of pheotypes and covariates) on a given cohort. You just need to provide:

- a models tables describing the models of interest
- a tables containing all the traits (pheotypes and covarites)
- a threshold for the maximum fraction of missing values in phenotype variables

Based on this information we first perform a preparing step that configure uniform analysis runs by grouping together compatible models. The process works as follows:

1. Phenotype values with the same data type and the same set of cvariates are grouped together
2. Within this groups pheotypes are further grouped to create multi-phenotype analysis in which none of the tested phenotype has more then the allowed fraction of missing values
3. The resulting analysis chunks are then automatically submitted to the pipeline for association analysis

Results for each analysis chunk are saved in a separate folder containing the configuration and the exact tables of phenotypes and covariates used in the execution.

In this way multiple hypotheses can be easily tested on large cohort with minimum amount of data pre-processing.

Note that in the output all phenotypes are renamed according to the model id defined in the input models table. In this way, each phenotype code in the output represents a unique model (combination of a phenotype and covariates). So for example, suppose you configured a multi-model run to test the phenotype `height` with 2 different sets of covariates:

| model_id | model |
| -------- | ----- |
| M1 | height ~ cov1+cov2 |
| M2 | height ~ cov3+cov4 |

The resulting output traits will be named `M1` and `M2` and the output files containing summary statistics will be named `M1.regenie.gz`, `M2.regenie.gz`.

See the [multi models configuration section](input-models-table.md) for more details on how to prepare the input files for a multi-models run.
