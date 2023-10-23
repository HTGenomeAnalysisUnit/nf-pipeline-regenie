# Models table

When executing a multi-models run you need to provide all phenotypes and covariates in a single file using the `phenotypes_filename` parameter. The exact data types and column names for phenotypes and covariates are then configured in the models table. 

The parameter `models_table` is used to activate a multi-models run by configuring a models table. The models table is a tab-separated file with header and columns names: model_id, model, trait_type, genetic_model, cat_var. 

- model_id: a unique identifier for the model.
- model: model descrition using col names from the phenotype table in the form `phenotype ~ covar1 + covar2 + ...`. Use `phenotype ~ 1` if you don't have any covariate
- trait_type: 'log' (for binary phenotype) or 'quant' (for quantitative phenotype).
- genetic_model: 'additive', 'dominant' or 'recessive'.
- cat_var: comma-separated list of categorical covariates. Use NA if none is present. Note that a max of 10 levels are accepted for a catagorical covariates by default, but this can be changed using the `maxCatLevels` parameter.

An example of a models table is shown below:

```bash
model_id        model   trait_type      genetic_model cat_var
M1      QP1 ~ Q1+Q2+Q5  quant   additive  NA
M2      QP2 ~ Q1+Q2+Q5  quant   additive  Q1,Q2
M3      QP3 ~ Q1+Q2+Q5  quant   additive  NA
M4      QP4 ~ 1  quant   additive   NA
```
