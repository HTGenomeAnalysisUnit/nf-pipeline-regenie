# nf-pipeline-regenie

A nextflow pipeline to perform genome-wide association studies (GWAS) using [regenie](https://github.com/rgcgithub/regenie) at high speed.

Original concept based on this amazing [github repository](https://github.com/genepi/nf-gwas) from Institute of Genetic Epidemiology, Innsbruck maintained by Sebastian Schönherr and Lukas Forer.

The pipeline is optimized for massive scaling and can analyze 5 quantitative phenotype on 500k individuals and 11M SNPs in ~4.5h.

## Documentation
Documentation can be found --DOCS IN PROGRESS--.

## Quick Start

1) Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (>=21.04.0)

2) Run the pipeline on a test dataset

```
nextflow run https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie -r v0.1 \
    -profile test,<slurm>
```

3) Run the pipeline on your data

```
nextflow run https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie -r v0.1 \
    -profile slurm -c your_project.conf
```

Please click [here](tests) for available test config files.

## Set up your project

An example of project config file is provided as `example_project.conf`. You can copy this file and adapt it to reflect your project data.

Essential inputs are:
1. A unique project id (`project`) that will be used to create output folder and reports
2. A bgen file of your full genotype data (`genotypes_imputed`). Regenie step2 will run faster on bgen v1.2 with 8 bits encoding. You can convert existing data using plink2 with `--export bgen-1.2 'bits=8'` option. No QC is performed on this file.
3. A bgi index for you bgen file. For a dataset named `my_dataset.bgen` the expected name of index is `my_dataset.bgen.bgi`. You can generate this using [bgenix tool](https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk). A `bgen` module is available to index your file if needed. You can use something like `bgenix -g my_data.bgen -index`
4. A SNP list for your full genotype data. This is a tab-separated file without header containing 2 columns: chr, pos for all SNPs in your data. You can obtain it from bgi index or from a bim file 
5. A {bed,bim,fam} dataset containing independent SNPs used by regenie step 1 (`genotypes_array`). Ideally, this would contain ~500k QCed SNPs and MUST contain less than 1M SNPs. An automatic QC will be performed on this file since step1 requires strict filtering criteria.
6. A tab-separated file with header for phenotypes and the list of column names to process. Note that quantitative and binary traits can not be mixed. Regenie will impute missing values when present.
7. A tab-separated file with header for covariates and the list of column names to process(can be omitted if no covars). Here it is fine to mix binary and quantitative covars. Binary covars must contain 0/1 values to be treated as binary. Note that **no missing values** are allowed in this table.

**NB.** The pipeline can also accept `vcf` file as input for full genotype data, and can generate `bgi` index and `snplist` file if missing. Note that in case of a large dataset, this will add considerable time to the execution due to slow conversion so it is strongly suggested to pre-process your input dataset to generate the needed inputs (BGEN + BGI + SNPLIST).

## Ouputs

By default the pipeline will generate all results in `output/<project>` organized as follow:
```
|-- logs
|-- regenie_step1_preds
|   |-- regenie_step1_out_1.loco.gz
|   `-- [one per test phenotype]
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

## License
nf-highspeed-gwas is MIT Licensed.

## Contact
* [Edoardo Giacopuzzi](mailto:edoardo.giacopuzzi@fht.org)
