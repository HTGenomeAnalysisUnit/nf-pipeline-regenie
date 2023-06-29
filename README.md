# nf-pipeline-regenie

A nextflow pipeline to perform genome-wide association studies (GWAS) and rare variant association analysis using [regenie](https://github.com/rgcgithub/regenie) at high speed.

Original concept based on this amazing [github repository](https://github.com/genepi/nf-gwas) from Institute of Genetic Epidemiology, Innsbruck maintained by Sebastian Schönherr and Lukas Forer.

The pipeline is optimized for massive scaling. In our tests it can analyze 100 quantitative phenotypes on teh full UKBB dataset (about 500k individuals) and 40M SNPs in ~4h under normal HPC load considering a limit of 250 concurrent tasks. When computational resources are available, you can cut down run time by increasing the limit on concurrent tasks.

Setting `genotypes_imputed` input trigger the GWAS analsysi, while `genotypes_rarevar` trigger the rare variant analysis using burden test and the set of test configured `rarevars_vc_test`. Possible tests include: skat, skato, acato, acatv, skato-acat.

Two running modes are available: **single project mode** and **multi models mode**.

## Table of contents

- [nf-pipeline-regenie](#nf-pipeline-regenie)
  - [Table of contents](#table-of-contents)
  - [Quick Start](#quick-start)
  - [Prepare config files](#prepare-config-files)
    - [Main parameters to adjust](#main-parameters-to-adjust)
  - [Understand input files for genetic data](#understand-input-files-for-genetic-data)
    - [Full genotype data (from imputation or sequencing) - MANDATORY](#full-genotype-data-from-imputation-or-sequencing---mandatory)
      - [vcf format](#vcf-format)
      - [bgen format](#bgen-format)
      - [pgen format](#pgen-format)
      - [bed format](#bed-format)
      - [Input dataset split by chromosome](#input-dataset-split-by-chromosome)
    - [QCed genotyped SNPs - MANDATORY](#qced-genotyped-snps---mandatory)
    - [Variant annotation files - MANDATORY for RARE VARIANT ANALYSIS](#variant-annotation-files---mandatory-for-rare-variant-analysis)
    - [LD panel - OPTIONAL (recommended for very large datasets)](#ld-panel---optional-recommended-for-very-large-datasets)
  - [Run in single project mode](#run-in-single-project-mode)
    - [Inputs for single project mode](#inputs-for-single-project-mode)
  - [Run in multi models mode](#run-in-multi-models-mode)
    - [Inputs for multi models mode](#inputs-for-multi-models-mode)
    - [Multi models execution monitoring](#multi-models-execution-monitoring)
    - [Resume execution for a failed task](#resume-execution-for-a-failed-task)
    - [Note on unexpected exit from multi models execution](#note-on-unexpected-exit-from-multi-models-execution)
  - [Important notes - Please read this](#important-notes---please-read-this)
    - [Phenotypes](#phenotypes)
    - [MAC filtering](#mac-filtering)
    - [Clumping](#clumping)
  - [Monitor execution on Nextflow tower](#monitor-execution-on-nextflow-tower)
  - [Annotation of tophits](#annotation-of-tophits)
  - [Outputs](#outputs)
  - [Re-use Step 1 predictions](#re-use-step-1-predictions)
  - [Full parameters explanation](#full-parameters-explanation)
  - [License](#license)
  - [Contact](#contact)


## Quick Start

1. Create a folder for your project (e.g. `yourproject`) 

2. Prepare the required genetic data for step 2, usually and [imputed dataset](#full-genotype-data-from-imputation---mandatory), and step 1, usually a [QCed genotyped dataset](#qced-genotyped-snps---mandatory). Then see the instruction to prepare config files for [single project run](#run-in-single-project-mode) or [multi models run](#run-in-multi-models-mode).

3. Prepare the [config files](#prepare-config-files) for your project in your project folder.

4. Invoke the pipeline

   Usually, you want to prepare a script to submit the pipeline in your project folder. In this example we use `sbatch` submission system, but this can be adapted to any scheduler:

   ```bash
   #!/bin/bash
   #SBATCH --job-name nf-regenie
   #SBATCH --output nf-regenie_master_%A.log
   #SBATCH --partition cpuq
   #SBATCH --cpus-per-task 1
   #SBATCH --mem 4G
   #SBATCH --time 1-00:00:00

   module load nextflow/22.10.1 singularity/3.8.5

   ###  For a single project  ###
   export NXF_OPTS="-Xms1G -Xmx4G" 
   nextflow run HTGenomeAnalysisUnit/nf-pipeline-regenie \
      -profile singularity,ht_cluster -c your_project.conf

   ###  For multiple models using a model table  ###
   nextflow run HTGenomeAnalysisUnit/nf-pipeline-regenie \
      -profile singularity,ht_cluster -c shared_parameters.conf \
      --with_master \
      --shared_config_file shared_parameters.conf \
      --models_table models.tsv \
      --traits_table full_traits.csv \
      --master_outdir outputs
   ```


Optionally, you can clone the latest pipeline version into the project folder using

`git clone --depth 1 https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie.git`

This will create a new folder called `nf-pipeline-regenie` in the current folder containing all the pipeline files.

You can eventually chose a specific version of the pipeline using the `--branch` option

`git clone --depth 1 --branch v1.7.4 https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie.git`

## Prepare config files

To run the pipeline, you need to prepare a config file. The following config file can be used as example:

- To run a single project, copy the template file `templates/single_project.conf` and adapt the parameters according to your input files and preferences.
- To run multiple models, copy the template file `templates/shared_parameters.conf` and adapt the parameters according to your input files and preferences.
- To create a profile for your HPC cluster, copy the template file `templates/profile_template.config` and adapt the parameters according to your HPC cluster configuration. Then use the profile in the pipeline execution by providing the new config file and the profile name, something like `-c myprofile.config -profile <profile_name>`.
  
### Main parameters to adjust

- If you provide a project id (`project`) that will be used in reports and a folder will be created with the same name to store all results.

- Set `chromosomes` to represent the list of chromosome to be included in the analysis

- Set `genotypes_build` to the build of your genotype data, either hg19 or hg38.

- Set `step2_gwas_chunk_size` and `step2_rarevar_chunk_size` according to the size of your dataset. These control how the dataset is split for step2 analysis. Teh default values usually work fine, but you can increase/decrease them if you have very large or small datasets. Keep in mind that a small chunk size may result in a very large amount of parallel jobs. By default the pipeline job submission rate is limited to 250 concurrent jobs. The total number of jobs will be *N_SNPs/gwas_chunk_size* for GWAS and *N_genes/rare_chunk_size* for rare variant analysis.

- Set `regenie_test` to 'additive', 'dominant' or 'recessive' according to the genetic model you want to test.

- Set `regenie_gwas_min_mac` and `regenie_rarevar_min_mac` to control the min allowed MAC for variants tested in either GWAS or rare variant analysis. Variants with MAC below this threshold will be excluded from the analysis.

- Set `annotation_min_log10p` to the min value ( -log10(pval) ) for top hit SNPs from GWAS results. These SNPs are also annotated in the manhattan plot in the HTML report.

- Set `rarevar_min_log10p` to the min value ( -log10(pval) ) for top hit genes from rare variants analysis.           

- Set `rarevar_tophits_min_value` and `rarevar_stat_test` to customise genes annotated in the manhattan plot in HTML report. For `rarevar_stat_test` you can select "pvalue", "FDR_bygroup", "FDR_alltests", "BONF_bygroup", "BONF_alltests". Then genes with ( -log10(stat_test) ) > min_value will be highlighted in the HTML report.

- Set `clump_p1` to the maximum pvalue allowed for index SNPs during plink clumping to define top loci

- If you are analyzing many phenotypes, it may be useful to set `make_report` to false. Generating the HTML graphical report can add a considerable amount of time for large datasets with many millions SNPs thus slowing down the overall execution.

- If you have categorical covariates, the maximum number of allowed categories is set to 10 bu default. You can adjsut this using the `maxCatLevels` parameter.

## Understand input files for genetic data

### Full genotype data (from imputation or sequencing) - MANDATORY

Regenie step2 will run faster on bgen v1.2 with 8 bits encoding. You can convert existing data using plink2 with `--export bgen-1.2 'bits=8'` option. No QC is performed on input data so ensure it is clean.

Input dataset for GWAS analysis is defined by `genotypes_imputed` and `genotypes_imputed_format` parameters, while input dataset for rare variants analysis is defined by `genotypes_rarevar` and `genotypes_rarevar_format`. The format parameter accepts 4 possible formats: vcf, bgen, pgen and bed.

The input dataset can be provided as a single file or split by chromosome including a `{CHROM}` tag in the input filename (see below).

Note that the pipeline will do highly parallelized access to the input dataset. Thus, especially if you are using a single file as input and when step 2 chunk size is small, it is reccomended to use a high performance filesystem optimized for parallel access.

#### vcf format

- input format: `'vcf'`
- input dataset: `'path/to/my_data.vcf.gz'`

The input VCF will be converted to `pgen` format using plink2 and `--double-id` option. The following parameters are used to control the conversion (see [relevant plink2 documentation](https://www.cog-genomics.org/plink/2.0/input#vcf)):

- `gwas_read_dosage_from` / `rarevar_read_dosage_from`: these params set from which field to read GT probabilities when converting VCF for the GWAS and rare variants input datasets, respectively. Accepted options are `'HDS'` (default) which usually works for VCF from Minimac4 imputation, `'DS'` for Minimac3 dosages or `'GP'` for genotype probabilities (to use with VCF from sequencing). Default is `'HDS'`.
- `import_dosage_certainty`: when `read_dosage_from = 'GP'` this parameter controls the minimum probability accepted to set a genotype. You can set to null to remove this filter. Default is `0.7`.
- `vcf_fixed_fid`: when `null` the conversion is performed using `--double-id` option, so that sample IDs in the VCF are used also for FID. If you set a string (like `'0'`) than the FID column in the VCF is fixe and set to this value for all samples. Default is `null`.

The converted dataset is saved to output folder when `save_pgen = true` (default).

**NB.** There are some aspects to keep in mind about the conversion from VCF to pgen format (see [plink2 docs](https://www.cog-genomics.org/plink/2.0/input#vcf) for more details):
- PGEN can store genotype dosage, but not genotype probabilities, thus probabilities are collapsed to dosages during conversio, which is usually fine for most analysis.
- When using `DS` or `HDS` during import, the dosage are read directly from the VCF.
- When using `GP`, the genotype probabilites are converted to dosages according to the `import_dosage_certainty` parameter. The default value of `0.7` means that a genotype dosage are set only if the probability of the most likely genotype is >= 0.7, otherwise teh whole genotype is set to missing.

#### bgen format

- input format: `'bgen'`
- input dataset: `'path/to/my_data.bgen'`

Some additional files are expected when using bgen format:

1. A bgi index for you bgen file. For a dataset named `my_dataset.bgen` the expected name of index is `my_dataset.bgen.bgi`. You can generate this using [bgenix tool](https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk). You can use something like `bgenix -g my_data.bgen -index`. 
2. A sample file. For a dataset named `my_dataset.bgen` the expected name of index is `my_dataset.sample`. This is standard sample file defined for the bgen format which contains sample level information. 
3. A SNP list for your dataset. This is a tab-separated file without header containing 6 columns: chr, id, cm, pos, ref, alt for all SNPs in your data, with chromosome in column 1 and position in column 4. Given a dataset named `my_dataset.bgen` the expected name of the snplist is `my_dataset.snplist`. You can generate a snplits from bgi index using

   ```bash
   bgenix -g my_dataset.begn -list | tail -n+3 \
   | awk '{OFS="\t"}; {print $3, $2, 0, $4, $6, $7}' | sed '$d' > my_dataset.snplist
   ```

**NB:** When using BGEN input, make sure that the sample ID in the BGEN or sample file can match FID + IID present in the covariates and phenotype input files, otherwise the pipeline will fail. Using a `.sample` file can help to have better control on the sample IDs. 

If any of these files is missing, the pipeline will generate them automatically and save them in the output folder by default. This behaviour can be controlled by `save_bgen_index`, `save_bgen_sample`, `save_snplist` parameters. Keep in mind that these steps can add a significant amount of time to the overall execution, so it is suggested to prepare these files in advance.

#### pgen format

- input format: `'pgen'`
- input dataset: `'path/to/my_data'`

For pgen input, you have to specify only the basename of the dataset. Given a input dataset named `my_dataset`, the pipeline will look for the following files: `my_data.pgen`, `my_data.pvar`, `my_data.psam`.

#### bed format

- input format: `'bed'`
- input dataset: `'path/to/my_data'`

For bed input, you have to specify only the basename of the dataset. Given a input dataset named `my_dataset`, the pipeline will look for the following files: `my_data.bed`, `my_data.bim`, `my_data.fam`.

#### Input dataset split by chromosome

If your input dataset is split by chromosome across multiple files, you can use the `{CHROM}` tag in your input file name. This tag must be placed corresponding to the number of chromosome in the filename. When using this method be careful that the chromosome names captured from the filename correspond to numbers 1-22 for autosomes. 

For example, if you have a dataset split by chromosome in the following way: 

- `my_dataset_chr1.bgen`
- `my_dataset_chr2.bgen`
- `my_dataset_chr3.bgen`
- ...

You can specify the input dataset as `my_dataset_chr{CHROM}.bgen` and the pipeline will automatically replace `{CHROM}` with the chromosome number.

### QCed genotyped SNPs - MANDATORY

A {bed,bim,fam} dataset containing independent SNPs used by regenie step 1 (`genotypes_array`). Ideally, this would contain ~500k QCed and pruned SNPs and MUST contain less than 1M SNPs. An additional QC will be automatically performed on this file since step1 requires strict filtering criteria. This need to be specified as prefix, so if you a dataset named `my_dataset.bed/bim/fam` you should specify `my_dataset` as input.

### Variant annotation files - MANDATORY for RARE VARIANT ANALYSIS

When running a rare variant analysis additional variant annotation files are needed and must be provided using the following parameters:

- `rarevar_set_list_file`: set list file as defined in regenie docs. Essentially a tab- or space-separated file with 4 columns: the set/gene name followed, a chromosome, physical position for the set/gene, a comma-separated list of variants included in the set/gene. This file is used to define the variant sets to be tested. The chromosome names must be numbers 1-22 for autosomes.

- `rarevar_anno_file`: variant annotation file for regenie. A tab- or space-separated file with 3 columns: variant name, the set/gene name, a single annotation category (for example missense, LoF, ...). Variants not in this file will be assigned to a default "NULL" category. A maximum of 63 annotation categories (+NULL category) is allowed.

- `rarevar_mask_file`: mask definition file for regenie. A tab- or space-separated file with 2 columns: a mask name followed by a comma-seperated list of categories included in the mask.


### LD panel - OPTIONAL (recommended for very large datasets)

LD panel files (`ld_panel` parameter). If you are analysing a large dataset with more than 50k samples, to speed up LD computation for clumping we suggest to prepare by chromosome bed/bim/fam files with the same variants present in the full genotype data (input bgen) but only a subset of samples. Then you can specify a pattern to these files like `/ld_panel/chr{CHROM}_ld`. The `{CHROM}` is automatically substituted with number 1-23 when the pipeline is running. This is only processed when `clumping` option is active. If the `ld_panel` parameter is not set and clumping is active the pipeline will use the full genotype data to estimate LD. Note that this will result in very long run time for huge datasets so providing an LD panel is highly reccomended when sample size is above 50k. You can generate LD files using plink2 and the `--keep` option to extract a subset of unrelated samples.

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

## Run in single project mode

In this mode you run a single GWAS model on the provided genetic data given a table of phenotypes and a table of covars.

### Inputs for single project mode

6. A tab-separated file with header for phenotypes (`phenotypes_filename`) and the list of column names to process (`phenotypes_columns`). Note that quantitative and binary traits can not be mixed. Regenie will impute missing values when present. Set `phenotypes_binary_trait` to true or false according to the type of phenotypes.
7. A tab-separated file with header for covariates (`covariates_filename`) and the list of column names (`covariates_columns`) to process(can be omitted if no covars). Here it is fine to mix binary and quantitative covars. Binary covars must contain 0/1 values to be treated as binary. Note that **no missing values** are allowed in this table.

**NB.** The first two columns of phenotype and covariate files must the named `FID` and `IID` and contain ids matching those present in the genotype files. Only samples found in all tables will be processed.

## Run in multi models mode

In this mode you can specifify a general trait table and a model table that describes the models you want to test (pheno ~ covars). Given a missingness threshold, the pipeline will automatically generate a list of of inputs and run the corresponding GWAS analyses.

### Inputs for multi models mode

6. A shared config file describing the input datasets and general config parameters (see `templates/shared_parameters.conf` for an example).
7. A tab-separated file with header (`full_traits.csv` in the example above) and first column named `IID` containing all traits (phenotypes and covariates) that are needed for the analysis.
8. A tab-separated file with header (`models.tsv` in the example above) and columns: model_id, model, trait_type, genetic_model, cat_var.
   - model_id: a unique identifier for the model.
   - model: model descrition using col names from `full_traits.csv` in the form `phenotype ~ covar1 + covar2 + ...`. Use `phenotype ~ 1` if you don't have any covariate
   - trait_type: 'log' (for binary phenotype) or 'quant' (for quantitative phenotype).
   - genetic_model: 'additive', 'dominant' or 'recessive'.
   - cat_var: comma-separated list of categorical covariates. Use NA if none is present. Note that a max of 10 levels are accepted for a catagorical covariates

```bash
model_id        model   trait_type      genetic_model cat_var
M1      QP1 ~ Q1+Q2+Q5  quant   additive  NA
M2      QP2 ~ Q1+Q2+Q5  quant   additive  Q1,Q2
M3      QP3 ~ Q1+Q2+Q5  quant   additive  NA
M4      QP4 ~ 1  quant   additive   NA
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

When clumping is active, the pipeline will save clumped data and clumps with genes annotation in the `toploci` folder. Note that if there are multiple identical SNP IDs clumping will fail. To avoid this you can for example modify your SNP ids to include ref/alt alleles as follows: `[SNPID]_[A0]_[A1]`.

**NB.** If you are using LD panel files as described in the [LD panel section](#ld-panel---optional-recommended-for-very-large-datasets), please ensure that SNP ids are concordant between bed/bim/fam files in the LD_panel and the input  dataset for GWAS analysis, otherwise SNPs that are not found will be dropped from clumping report.

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

## Annotation of tophits

For top hits SNPs / loci the pipeline performs automatic annotation with overlapping genes and nearby genes. The maximum distance between the variant / locus and the sorrounding annotated gene is set with the parameter `annotation_interval_kb` (25kb by default). Genes located +/- defined kb from the SNP / locus are annotated as nearby genes.

Pre-made BED files for gene annotation are available for GRCh37 and GRCh38. These files are generated from the [GENCODE annotation](https://www.gencodegenes.org/) version 39. When performing annotation one can chose to annotate for all genes or only protein coding genes by setting the `genes_group` parameter to `all` or `protein_coding` respectively. By default, the pipeline will annotate for protein coding genes.

If you want to use a different annotation file, you can provide custom gene definition with the `genes_bed` and `genes_ranges` parameters. In this case the pre-made annotation files will not be used. For this you need to provide:

- `genes_bed`: a regular BED file (chromosome, start, end; zero-based), with gene names in the 4th column. Note that a name must be present for all genes.
- `genes_ranges`: essentially the same as above but space separated.

## Outputs

By default the pipeline will generate all results in a folder named according to `project` parameter. When you provide the `outdir` parameter, the output folder will be created as `<outdir>/<project>`. Main structure of the output folder is as follows:

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
|-- results
|   |-- gwas
|   |   |-- tophits
|   |   `-- toploci
|   `-- rarevar
|       `-- tophits
|-- step2_chunks
`-- validated_input
```

- main results from GWAS and rare variants tests are in `results`. For GWAS, top associated SNPs annotated for gene(s) overlap and nearby gene(s) are saved in the `tophits` sub-folder, and top loci after clumping are saved in `toploci` sub folder. For rare variants, top associated genes are saved in `tophits` sub folder.
- step 1 predictions in `regenie_step1_preds` can be reused for further analyses on the same dataset as long as the input bgen, phenotype file and covars file are exactly the same and phenotypes and covars list are provided in the same order.
- logs from all operations are saved in `logs` for debugging. Note that more than a thousand log file may be generated when the input dataset is large.
- HTML reports for GWAS and rare variants are saved in `reports` folder
- when you are running in multi models mode, one folder will be created for each run_id under the `master_outdir` folder.

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
