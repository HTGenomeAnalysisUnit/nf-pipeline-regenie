# Outputs

By default, the pipeline will generate all results in a folder named according to the `project` parameter: `<project>_output`. Main results of variant and gene based analyses are then stored in this folder using the `project` parameter (i.e. `<project>_output/<project>`). When you provide the `outdir` parameter, this will be used as the main output folder, thus the association results will be placed in `<outdir>/<project>`.

When you are running in multi-model or multi-project mode you will end up having multiple sub-folders, one for each project id.

Main structure of the output folder is as follows:

```bash
.
|-- analysis_config
|-- step2_chunks
|   |-- step2_dataset_autosomes.GWAS-chunks.txt
|   `-- step2_dataset_autosomes.genes-chunks.txt
`-- project1
    |-- logs
    |-- regenie_step1_preds
    |-- reports
    |-- results
    |   |-- gwas
    |   `-- rare_var
    `-- validated_input
```

- validated input tables for phenotyopes and covariates are saved in `validated_input` folder. These are the actual tables used in the analysis after validation and formatting (sanitize column names, remove spaces, etc.).
- main results from GWAS and rare variants tests are in `results`. For GWAS, top associated SNPs annotated for gene(s) overlap and nearby gene(s) are saved in the `tophits` sub-folder, and top loci after clumping are saved in `toploci` sub folder. For rare variants, top associated genes are saved in `tophits` sub folder.
- step 1 predictions in `regenie_step1_preds` can be reused for further analyses on the same dataset as long as the input bgen, phenotype file and covars file are exactly the same and phenotypes and covars list are provided in the same order.
- general configuration files for the pipeline execution are saved in `analysis_config` folder.
- logs from all project-specific operations are saved in `logs` for debugging.
- HTML reports for GWAS and rare variants are saved in `reports` folder when `make_reports` is true.
