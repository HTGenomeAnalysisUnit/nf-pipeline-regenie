# Outputs

By default, the pipeline will generate all results in a folder named `nf-regenie-pipeline_results`. Main results of variant analysis are then stored in this folder using the `project` parameter (i.e. `<outdir>/<project>`). When you provide the `outdir` parameter, this will be used as the main output folder.

When you are running in multi models mode, one folder will be created for each run_id under the main output folder.

Main structure of the output folder is as follows:

```bash
.
|-- analysis_config
|-- logs
|   |-- clumping
|   |   `-- <pheno>_clump
|   |-- step2_gwas_logs
|   `-- step2_rarevar_logs
|-- regenie_step1_preds
|-- reports
|   |-- gwas
|   `-- rarevar
|-- results
|   |-- gwas
|   |   |-- tophits
|   |   `-- toploci
|   `-- rarevar
|       `-- tophits
|-- step2_chunks
`-- validated_input
```

- validated input tables for phenotyopes and covariates are saved in `validated_input` folder. These are the actual tables used in the analysis after validation and formatting (sanitize column names, remove spaces, etc.).
- main results from GWAS and rare variants tests are in `results`. For GWAS, top associated SNPs annotated for gene(s) overlap and nearby gene(s) are saved in the `tophits` sub-folder, and top loci after clumping are saved in `toploci` sub folder. For rare variants, top associated genes are saved in `tophits` sub folder.
- step 1 predictions in `regenie_step1_preds` can be reused for further analyses on the same dataset as long as the input bgen, phenotype file and covars file are exactly the same and phenotypes and covars list are provided in the same order.
- logs from all operations are saved in `logs` for debugging and the config files used for each run are saved in `analysis_config`.
- HTML reports for GWAS and rare variants are saved in `reports` folder.
