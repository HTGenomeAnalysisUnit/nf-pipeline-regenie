# Re-use step 1 predictions

If you set `save_step1_predictions` to true, the pipeline will save the step 1 predictions in `regenie_step1_preds` folder. 

These predictions can be re-used in  further analyses on the same dataset when phenotype(s) and covariates are the same. This is useful for example to perform a rare variant analysis on the same dataset used for a previous GWAS analysis.

You can load level 1 preds from this folder in subsequent analyses by setting `regenie_premade_predictions` to a path like `/results/regenie_step1_preds/regenie_step1_out*`.

The pipeline then expects the following:

- A file named `regenie_step1_out_pred.list` must be present
- One file per phenotype is expected named `regenie_step1_out_1.loco.gz` `regenie_step1_out_2.loco.gz`, ...
- phenotypes and covariates used in the new analyses must be exactly the same used to generate step1 predictions and column names for phenotypes must be configured in `phenotypes_columns` in the same order as in the original analysis.