# nf-pipeline-regenie

A nextflow pipeline to perform genome-wide association studies (GWAS) using [regenie](https://github.com/rgcgithub/regenie) at high speed.

Original concept based on this amazing [github repository](https://github.com/genepi/nf-gwas) from Institute of Genetic Epidemiology, Innsbruck maintained by Sebastian Schönherr and Lukas Forer.

The pipeline is optimized for massive scaling and can analyze 10 quantitative phenotype on 500k individuals and 40M SNPs in ~4.5h.

Two running modes are available: **single project mode** and **multi models mode**.

## Quick Start

1. Create a folder for your project (e.g. `yourproject`) and clone this pipeline version into the project folder using

`git clone --depth 1 --branch v1.0 https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie.git`

This will create a new folder called `nf-pipeline-regenie` in the current folder containing all the pipeline files.

2. Prepare the [required genetic data](#required-inputs) and the additional input tables needed for [single project](#single-project-mode) or [multi models](#multi-models-mode).

3. Prepare the [config files](#prepare-config-files) for your project.

4. Invoke the pipeline

   Ideally, you should prepare a script to submit the pipeline in your project folder using `sbatch`. The following template can be used as example:

   ```bash
   #!/bin/bash
   #SBATCH --job-name nf-regenie
   #SBATCH --output nf-regenie_master_%A.log
   #SBATCH --partition cpuq
   #SBATCH --cpus-per-task 1
   #SBATCH --time 15-00:00:00

   module load nextflow/21.10.6 singularity/3.6.3

   ###  For a single project  ###
   nextflow run nf-pipeline-regenie \
      -profile slurm -c your_project.conf

   ###  For multiple models using a model table  ###
   nextflow run nf-pipeline-regenie \
      -profile slurm -c shared_parameters.conf \
      --with_master \
      --shared_config_file shared_parameters.conf \
      --models_table models.tsv \
      --traits_table full_traits.csv \
      --master_outdir outputs
   ```

## Prepare config files

To run the pipeline, you need to prepare a config file. The following config file can be used as example:

- To run a single project, copy the template file `templates/single_project.conf` and adapt the parameters according to your input files and preferences.
- To run multiple models, copy the template file `templates/shared_parameters.conf` and adapt the parameters according to your input files and preferences.

### Main parameters to adjust

- If you provide a project id (`project`) that will be used in reports and a folder will be created with the same name to store all results.

- Set `step2_chunk_size` according to the size of your dataset. A value of 100000 usually works fine, but you can decrease this down to 25000 for very large datasets to speed up the process if you have large number of cpus available. Note that the pipeline will submit *N_SNPs/chunk_size* jobs in batches of 500 jobs.

- Set `regenie_test` to 'additive', 'dominant' or 'recessive' according to the genetic model you want to test.

## Required input files

### Full genotype data (from imputation)

1. A bgen file of your full genotype data (`genotypes_imputed`). Regenie step2 will run faster on bgen v1.2 with 8 bits encoding. You can convert existing data using plink2 with `--export bgen-1.2 'bits=8'` option. No QC is performed on this file so ensure it is clean.
2. A bgi index for you bgen file. For a dataset named `my_dataset.bgen` the expected name of index is `my_dataset.bgen.bgi`. You can generate this using [bgenix tool](https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk). A `bgen` module is available on our HPC to index your file if needed. You can use something like `bgenix -g my_data.bgen -index`
3. A SNP list for your full genotype data. This is a tab-separated file without header containing 2 columns: chr, pos for all SNPs in your data. You can obtain it from bgi index using 

   ```bash
   bgenix -g input.bgen -list \
   | tail -n+3 | cut -f3,4 | sed '$d' > ${prefix}.snplist
   ```

**Other input formats.** The pipeline can also accept `vcf` file as input for full genotype data, and can generate `bgi` index and `snplist` file if missing. Note that in case of a large dataset, this will add considerable time to the execution due to slow conversion so it is strongly suggested to pre-process your input dataset to generate the needed inputs (BGEN + BGI + SNPLIST).

### QCed genotyped SNPs

4. A {bed,bim,fam} dataset containing independent SNPs used by regenie step 1 (`genotypes_array`). Ideally, this would contain ~500k QCed and pruned SNPs and MUST contain less than 1M SNPs. An additional QC will be automatically performed on this file since step1 requires strict filtering criteria. This need to be specified as prefix, so if you a dataset named `my_dataset.bed/bim/fam` you should specify `my_dataset` as input.

### LD panel

5. LD panel files (`ld_panel` parameter). To speed up LD computation for clumping we suggest to prepare by chromosome bed/bim/fam files with the same variants present in the full genotype data (input bgen) but only a subset of samples. Then you can specify a pattern to these files like `/ld_panel/chr{CHROM}_ld`. The `{CHROM}` is automatically substituted with number 1-23 when the pipeline is running. If clumping is not active this can be omitted. If the `ld_panel` parameter is not set and clumping is active the pipeline will use the full genotype data to estimate LD. Note that this will result in very long run time for huge datasets so providing an LD panel is highly reccomended when sample size is above 20k. You can generate LD files from input BGEN using plink2

```bash
chrom=1 #Chromosome id
samples="samples.tsv" #File with a subset of unrelated individuals

plink2 \
--bgen input.bgen ref-first \
--sample input.sample \
--keep ${samples} \
--chr ${chrom} \
--make-bed \
--out chr${chrom}_ld_panel
```

## Single project mode

In this mode you run a single GWAS model on the provided genetic data given a table of phenotypes and a table of covars.

### Additional inputs for single project mode

6. A tab-separated file with header for phenotypes (`phenotypes_filename`) and the list of column names to process (`phenotypes_columns`). Note that quantitative and binary traits can not be mixed. Regenie will impute missing values when present. Set `phenotypes_binary_trait` to true or false according to the type of phenotypes.
7. A tab-separated file with header for covariates (`covariates_filename`) and the list of column names (`covariates_columns`) to process(can be omitted if no covars). Here it is fine to mix binary and quantitative covars. Binary covars must contain 0/1 values to be treated as binary. Note that **no missing values** are allowed in this table.

**NB.** The first two columns of phenotype and covariate files must the named `FID` and `IID` and contain ids matching those present in the genotype files. Only samples found in all tables will be processed.

## Multi models mode

In this mode you can specifify a general trait table and a model table that describes the models you want to test (pheno ~ covars). Given a missingness threshold, the pipeline will automatically generate a list of of inputs and run the corresponding GWAS analyses.

### Files for multi models mode

6. A shared config file describing the input datasets and general config parameters (see `templates/shared_parameters.conf` for an example).
7. A tab-separated file with header (`full_traits.csv` in the example above) and first column named `IID` containing all traits (phenotypes and covariates) that are needed for the analysis.
8. A tab-separated file with header (`models.tsv` in the example above) and columns: model_id, model, trait_type, genetic_model. 
   - model_id: a unique identifier for the model.
   - model: model descrition using col names from `full_traits.csv` in the form `phenotype ~ covar1 + covar2 + ...`.
   - trait_type: 'log' or 'quant'.
   - genetic_model: 'additive', 'dominant' or 'recessive'.

```bash
model_id        model   trait_type      genetic_model
M1      QP1 ~ Q1+Q2+Q5  quant   additive
M2      QP2 ~ Q1+Q2+Q5  quant   additive
M3      QP3 ~ Q1+Q2+Q5  quant   additive
M4      QP4 ~ Q1+Q2+Q5  quant   additive
```

### Multi models execution monitoring

If you monitor the execution using nextflow, the master job will be named `nf-highspeed-gwas` and the single GWAS jobs will be named `run-<chunk_id>-gwas`.

In the master output folder you will also see:

- `master_config` folder containing logs, info and config about the master job
- `execution_log` folder. This will contain the exection path and status for each job and eventually error log files.
  
### Resume execution for a failed task

When one of the task fails, the pipeline will report a warning, but continue to execute the remaining tasks. When the pipeline is finished, you can inspect `job_execution_summary.log` located in `<output_dir>/execution_log` and check which task have failed and the corresponding running directory. An `_error.log` file will be present in the same directory explainign the error.

After fixing the error, you can resume execution of a specific task following these steps:

1. move into the corresponding running directory and add the `-resume` flag to the command in the `.command.sh` file.
2. remove the setting folder from the directory (it is named `chunk_XX`)
3. submit your job again by simply running `sbatch .command.run`.

### Note on unexpected exit from multi models execution

At the moment, if the master pipeline terminates unexpectedly, it is likely that jobs realted to the single model executions will not be cleaned up. This is because the master pipeline is not aware of the jobs related to the single model executions.
In this case, please verify if you have any running jobs using `squeue -u $USER` and terminate any job with name containing `nf-SETUP_MULTIPLE_RUNS_SUBMIT_GWAS_RUN` or `nf-NF_GWAS_REGENIE`.

Normally, if one of the single run submission terminates with error, the master pipeline will go on and a warning is reported. This ensure the whole pipeline can gracefully terminate and all sub-jobs are properly cleaned up. You can see from the warning message which run has failed an eventually re-run this analysis as single GWAS.

## Important notes - Please read this

### Phenotypes

- Phenotype file must have header and first 2 columns must be FID, IID
- Samples listed in pheno file that are not in bgen/bed/pgen file are ignored. Genotyped samples that are not in the pheno file are removed from the analysis.
- With quantitative traits, missing values are mean-imputed in Step 1 and they are dropped when testing each phenotype in Step 2
- With binary traits, missing values are mean-imputed in Step 1 when fitting the level 0 linear ridge regression and they are dropped when fitting the level 1 logistic ridge regression for each trait. In Step 2, missing values are dropped when testing each trait.

### MAC filtering

Note that the pipeline impose a MAC / MAF filter on genotyped data used for step1, as reccomended by regenie developers, and a min MAC 50 at step2 as suggested in the recent large proteome analysis of UKBB (see [the preprint](https://www.biorxiv.org/content/10.1101/2022.06.17.496443v1)).

This MAC filters are applied after subsetting the data to contain only samples seen in the phenotype table. FID/IID are extracted from the pheno file and used with `--keep` in genotyped data QC to generate a QCed dataset used in regenie step 1. Then in step 2, regenie automatically removes any sample not seen in step 1 predictions before applying the MAC filter.

As a result of this process, if you use only a subset of samples present in the original imputed data, it is possible that some variants will be filtered out due to low MAC. Check log files and the results table to see exactly how many variants were considered in the analysis.

### Clumping

When clumping is active, the pipeline will save clumped data and clumps with genes annotation in the `toploci` folder. Note that SNP IDs in these files are made from original IDs as follows: `[SNPID]_[A0]_[A1]` to avoid duplicate conflict when there are multiple SNPs with the same rsID, but different tested alleles.

**NB.** Please ensure that alleles order is concordant between bed/bim/fam files in the LD_panel and the BGEN file of imputed data. See [LD panel](#ld-panel) section.

## Monitor execution on Nextflow tower

To easily monitor pipeline execution, we suggest to use Nextflow tower. First, register to the [Nextflow tower](https://cloud.tower.nf/) using your GitHub or Google account. Then, click on your profile (upper right corner) and select `Your tokens`. Then click add token and follow the instructions to create a new token. Make sure to copy the generated token, since you were not able to see it again.

Finally, add the following to `shared_parameters.conf` or `single_project.conf` file:

```nextflow
tower {
  enabled = true
  accessToken = 'your_token_here'
}
```

Now when the pipeline is running you should be able to monitor progress in the [Nextflow tower](https://cloud.tower.nf/) under your runs.

## Outputs

By default the pipeline will generate all results in a folder named according to `project` variable. Alternatively you can provide an `outdir` parameter and the output folder will be created as `<outdir>/<project>`. Output files are organized as follows:

```bash
|-- logs
|-- analysis_config
|-- regenie_step1_preds
|   |-- regenie_step1_out_1.loco.gz
|   |-- [one per tested phenotype]
|   `-- regenie_step1_out_pred.list
|-- results  
|   |-- QP1.regenie.gz
|   |-- QP1.regenie.gz.tbi
|   |-- [one pair per test phenotype] 
|   |-- tophits
|   `-- toploci [only when clumping is active]
|-- step2_chunks
|   `-- chunks.txt
|-- test-gwas-split.QP1.regenie.html
|-- [one html report per tested phenotype]
`-- validated_input
    |-- covars.validated.txt
    `-- phenos.validated.txt
```

- main results from association are in `results`, with top associated SNPs annotated for gene(s) overlap and nearby gene(s) located in the `tophits` sub-folder
- step 1 predictions in `regenie_step1_preds` can be reused for further analyses on the same dataset as long as the input bgen, phenotype file and covars file are exactly the same and phenotypes and covars list are provided in the same order.
- logs from all operations are saved in `logs` for debugging. Note that more than a thousand log file may be generated when the input dataset is large.
- when you are running in multi models mode, one folder will be created for each run_id under the `master_outdir` folder.

## DB function

The pipeline include steps to generate a `bcf` based DB storing association results and models details. With default only results with P < 0.05 are stored, but you can adjust this by setting `db_pval_threshold` to the desired -LOG10P value.

Once created DB from multiple runs can be merged easily and DB in this format allows fast access across results. In our test on ~10M SNPs, it can recover association values (P, EFFECT, SE) for a SNPs across 5000 traits (pheWAS) in few secs and recover all results in a 500kb window in 20-30 secs. The implementation is composed by a BCF file storing the actual association results and a SQLITE db file storing the SNPs informations (coordinate, unique marker ID and rsID).

However, generating the DB can add significant time to the overall execution so this function is disabled by default. You can activate it by setting `create_db` to true.

## Re-use Step 1 predictions

If you set `save_step1_predictions` to true, the pipeline will save the step 1 predictions in `regenie_step1_preds` folder. This can be used to re-use the predictions for further analysis.
You can load level 1 preds from this folder in subsequent analyses by setting `regenie_premade_predictions` to a path like `/results/regenie_step1_preds/regenie_step1_out*`.

- A file named regenie_step1_out_pred.list must be present
- One file per phenotype is expected named `regenie_step1_out_1.loco.gz` `regenie_step1_out_2.loco.gz`, ...
- phenotypes and covariates used in the new analyses must be exactly the same used to generate step1 predictions (both files and column designation must match exactly)

## Full parameters explanation

See the [project wiki](https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie/-/wikis/home) for detailed explanation of all configuration parameters

## License

nf-highspeed-gwas is MIT Licensed.

## Contact

- [Edoardo Giacopuzzi](mailto:edoardo.giacopuzzi@fht.org)
