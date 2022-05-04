# nf-pipeline-regenie

A nextflow pipeline to perform genome-wide association studies (GWAS) using [regenie](https://github.com/rgcgithub/regenie) at high speed.

Original concept based on this amazing [github repository](https://github.com/genepi/nf-gwas) from Institute of Genetic Epidemiology, Innsbruck maintained by Sebastian Schönherr and Lukas Forer.

The pipeline is optimized for massive scaling and can analyze 5 quantitative phenotype on 500k individuals and 11M SNPs in ~4.5h.

## Quick Start

1) Load Nextflow and Singularity modules
   
```bash
module load nextflow/22.04.0 singularity/3.6.3
```

2A) Run the pipeline on a single project

```bash
nextflow run /project/alfredo/pipelines/nf-pipeline-regenie/main.nf \
    -profile slurm -c your_project.conf
```

2B) Run the pipeline for multiple models using a model table

```bash
nextflow run /project/alfredo/pipelines/nf-pipeline-regenie/main.nf \
    -profile slurm -c shared_parameters.conf \
    --with_master \
    --shared_config_file shared_parameters.conf \
    --models_table models.tsv \
    --traits_table full_traits.csv \
    --master_outdir outputs
```

Ideally, you should prepare a script to submit the pipeline using `sbatch`. The following template can be used as example:

```bash
#!/bin/bash
#SBATCH --job-name nf-regenie
#SBATCH --output nf-regenie_master_%A.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1

module load nextflow/22.04.0 singularity/3.6.3

###  For a single project  ###
nextflow run /project/alfredo/pipelines/nf-pipeline-regenie/main.nf \
   -profile slurm -c your_project.conf

###  For multiple models using a model table  ###
nextflow run /project/alfredo/pipelines/nf-pipeline-regenie/main.nf \
   -profile slurm -c shared_parameters.conf \
   --with_master \
   --shared_config_file shared_parameters.conf \
   --models_table models.tsv \
   --traits_table full_traits.csv \
   --master_outdir outputs
```

## Set up configuration file

- To run a single project, copy the template file `templates/single_project.conf` and adapt the parameters.
- To run multiple models, copy the template file `templates/shared_parameters.conf` and adapt the parameters.

## Understand the inputs

### Shared inputs

1. A unique project id (`project`) that will be used to create output folder and reports
2. A bgen file of your full genotype data (`genotypes_imputed`). Regenie step2 will run faster on bgen v1.2 with 8 bits encoding. You can convert existing data using plink2 with `--export bgen-1.2 'bits=8'` option. No QC is performed on this file so ensure it is clean.
3. A bgi index for you bgen file. For a dataset named `my_dataset.bgen` the expected name of index is `my_dataset.bgen.bgi`. You can generate this using [bgenix tool](https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk). A `bgen` module is available to index your file if needed. You can use something like `bgenix -g my_data.bgen -index`
4. A SNP list for your full genotype data. This is a tab-separated file without header containing 2 columns: chr, pos for all SNPs in your data. You can obtain it from bgi index or from a bim file 
5. A {bed,bim,fam} dataset containing independent SNPs used by regenie step 1 (`genotypes_array`). Ideally, this would contain ~500k QCed and pruned SNPs and MUST contain less than 1M SNPs. An automatic QC will be performed on this file since step1 requires strict filtering criteria. This need to be specified as prefix, so if you a dataset named `my_dataset.bed/bim/fam` you should specify `my_dataset` as input.
6. Set `step2_chunk_size` according to the size of your dataset. A value of 25000 usually works fine, but you can decrease this down to 10000 for large datasets (> 300k individuals) or increase it up to 100000 for smaller ones (less than 50k individuals)

**Other input formats.** The pipeline can also accept `vcf` file as input for full genotype data, and can generate `bgi` index and `snplist` file if missing. Note that in case of a large dataset, this will add considerable time to the execution due to slow conversion so it is strongly suggested to pre-process your input dataset to generate the needed inputs (BGEN + BGI + SNPLIST).

### Files for single project mode

7. A tab-separated file with header for phenotypes (`phenotypes_filename`) and the list of column names to process (`phenotypes_columns`). Note that quantitative and binary traits can not be mixed. Regenie will impute missing values when present.
8. Set `phenotypes_binary_trait` to true or false according to the type of phenotypes.
9. A tab-separated file with header for covariates (`covariates_filename`) and the list of column names (`covariates_columns`) to process(can be omitted if no covars). Here it is fine to mix binary and quantitative covars. Binary covars must contain 0/1 values to be treated as binary. Note that **no missing values** are allowed in this table.
10. Set `regenie_test` to 'additive', 'dominant' or 'recessive' according to the genetic model you want to test.

**NB.** The first two columns of phenotype and covariate files must the named `FID` and `IID` and contain ids matching those present in the genotype files. Only samples found in all tables will be processed.

### Files for multi models mode

7. A shared config file describing the input datasets and general config parameters (see `templates/shared_parameters.conf` for an example).
8. A tab-separated file with header (`full_traits.csv` in the example above) and first column named `IID` containing all traits (phenotypes and covariates) that are needed for the analysis.
9. A tab-separated file with header (`models.tsv` in the example above) and columns: model_id, model, trait_type, genetic_model. 
   - model_id: a unique identifier for the model.
   - model: model descrition using col names from `full_traits.csv` in the form `phenotype ~ covar1 + covar2 + ...`.
   - trait_type: 'log' or 'quant'.
   - genetic_model: 'additive', 'dominant' or 'recessive'.

```
model_id        model   trait_type      genetic_model
M1      QP1 ~ Q1+Q2+Q5  quant   additive
M2      QP2 ~ Q1+Q2+Q5  quant   additive
M3      QP3 ~ Q1+Q2+Q5  quant   additive
M4      QP4 ~ Q1+Q2+Q5  quant   additive
```

## Outputs

By default the pipeline will generate all results in a folder named according to `project` variable. Alternatively you can provide an `outdir` parameter and the output folder will be created as `<outdir>/<project>`. Output files are organized as follows:
```
|-- logs
|-- regenie_step1_preds
|   |-- regenie_step1_out_1.loco.gz
|   |-- [one per tested phenotype]
|   `-- regenie_step1_out_pred.list
|-- results  
|   |-- QP1.regenie.gz
|   |-- QP1.regenie.gz.tbi
|   |-- [one pair per test phenotype] 
|   `-- tophits
|-- step2_chunks
|   `-- chunks.txt
|-- test-gwas-split.QP1.regenie.html
|-- [one html report per test phenotype]
`-- validated_input
    |-- covars.validated.txt
    `-- phenos.validated.txt
```

- main results from association are in `results`, with top associated SNPs annotated for gene(s) overlap and nearby gene(s) located in the `tophits` sub-folder
- step 1 predictions in `regenie_step1_preds` can be reused for further analyses on the same dataset as long as the input bgen, phenotype file and covars file are exactly the same and phenotypes and covars list are provided in the same order.
- logs from all operations are saved in `logs` for debugging. Note that more than a thousand log file may be generated when the input dataset is large.
- when you are running in multi models mode, one folder will be created for each run_id under the `master_outdir` folder.

## Monitor execution on Nextflow tower

To easily monitor pipeline execution, we suggest to use Nextflow tower. First, register to the [Nextflow tower](https://cloud.tower.nf/) using your GitHub or Google account. Then, click on your profile (upper right corner) and select `Your tokens`. Then click add token and follow the instructions to create a new token. Make sure to copy the generated token, since you were not able to see it again.

Finally, add the following to `shared_parameters.conf` or `single_project.conf` file:

```
tower {
  enabled = true
  accessToken = 'your_token_here'
}
```

Now when the pipeline is running you should be able to monitor progress in the [Nextflow tower](https://cloud.tower.nf/) under your runs.

## Re-use Step 1 predictions

If you set `save_step1_predictions` to true, the pipeline will save the step 1 predictions in `regenie_step1_preds` folder. This can be used to re-use the predictions for further analysis. 
You can load level 1 preds from this folder in subsequent analyses by setting `regenie_premade_predictions` to a path like `/results/regenie_step1_preds/regenie_step1_out*`.

- A file named regenie_step1_out_pred.list must be present 
- One file per phenotype is expected named `regenie_step1_out_1.loco.gz` `regenie_step1_out_2.loco.gz`, ...
- phenotypes and covariates used in the new analyses must be exactly the same used to generate step1 predictions (both files and column designation must match exactly)

## Important notes on handling phenos

- Samples listed in pheno file that are not in bgen/bed/pgen file are ignored. Genotyped samples that are not in this file are removed from the analysis.
- With quantitative traits, missing values are mean-imputed in Step 1 and they are dropped when testing each phenotype in Step 2
- With binary traits, missing values are mean-imputed in Step 1 when fitting the level 0 linear ridge regression and they are dropped when fitting the level 1 logistic ridge regression for each trait. In Step 2, missing values are dropped when testing each trait.

## Note on unexpected exit from multi models execution

At the moment, if the master pipeline terminates unexpectedly, it is likely that jobs realted to the single model executions will not be cleaned up. This is because the master pipeline is not aware of the jobs related to the single model executions.
In this case, please verify if you have any running jobs using `squeue -u $USER` and terminate any job with name containing `nf-SETUP_MULTIPLE_RUNS_SUBMIT_GWAS_RUN` or `nf-NF_GWAS_REGENIE`.

Normally, if one of the single run submission terminates with error, the master pipeline will go on and a warning is reported. This ensure the whole pipeline can gracefully terminate and all sub-jobs are properly cleaned up. You can see from the warning message which run has failed an eventually re-run this analysis as single GWAS.

## Documentation
More detailed documentation of all parameters can be found --DOCS IN PROGRESS--.

## License
nf-highspeed-gwas is MIT Licensed.

## Contact
* [Edoardo Giacopuzzi](mailto:edoardo.giacopuzzi@fht.org)
