# Full genotype data 

To perform genetic association or rare variants analysis (regenie step 2) you must provide the full genetic dataset for your study, usually a large genetic dataset from imputation or sequencing. No QC is performed on teh full genetic dataset so ensure it is clean.

Input dataset for GWAS analysis is defined by `genotypes_imputed` and `genotypes_imputed_format` parameters, while input dataset for rare variants analysis is defined by `genotypes_rarevar` and `genotypes_rarevar_format`. The format parameter accepts 4 possible formats: vcf, bgen, pgen or bed.

The input dataset can be provided as a single file or split by chromosome including a `{CHROM}` tag in the input filename (see below).

## Input formats

### vcf format

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

### bgen format

- input format: `'bgen'`
- input dataset: `'path/to/my_data.bgen'`

Some additional files are expected when using bgen format:

1. A bgi index for you bgen file. For a dataset named `my_dataset.bgen` the expected name of index is `my_dataset.bgen.bgi`. You can generate this using [bgenix tool](https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk). You can use something like `bgenix -g my_data.bgen -index`. 
2. A sample file. For a dataset named `my_dataset.bgen` the expected name of index is `my_dataset.sample`. This is standard sample file defined for the bgen format which contains sample level information. 
3. A SNP list for your dataset. This is a tab-separated file without header containing 6 columns: chr, id, cm, pos, ref, alt for all SNPs in your data, with chromosome in column 1 and position in column 4. Given a dataset named `my_dataset.bgen` the expected name of the snplist is `my_dataset.snplist`. You can generate a snplits from bgi index using

   ```bash
   bgenix -g my_dataset.bgen -list | tail -n+3 \
   | awk '{OFS="\t"}; {print $3, $2, 0, $4, $6, $7}' | sed '$d' > my_dataset.snplist
   ```

**NB:** When using BGEN input, make sure that the sample ID in the BGEN or sample file can match FID + IID present in the covariates and phenotype input files, otherwise the pipeline will fail. Using a `.sample` file can help to have better control on the sample IDs. 

If any of these files is missing, the pipeline will generate them automatically and save them in the output folder by default. This behaviour can be controlled by `save_bgen_index`, `save_bgen_sample`, `save_snplist` parameters. Keep in mind that these steps can add a significant amount of time to the overall execution for large datasets, so it is suggested to prepare these files in advance.

### pgen format

- input format: `'pgen'`
- input dataset: `'path/to/my_data'`

For pgen input, you have to specify only the basename of the dataset. Given a input dataset named `my_dataset`, the pipeline will look for the following files: `my_data.pgen`, `my_data.pvar`, `my_data.psam`.

### bed format

- input format: `'bed'`
- input dataset: `'path/to/my_data'`

For bed input, you have to specify only the basename of the dataset. Given a input dataset named `my_dataset`, the pipeline will look for the following files: `my_data.bed`, `my_data.bim`, `my_data.fam`.

## Input dataset split by chromosome

If your input dataset is split by chromosome across multiple files, you can use the `{CHROM}` tag in your input file name. This tag must be placed corresponding to the number of chromosome in the filename. When using this method be careful that the chromosome names captured from the filename correspond to numbers 1-22 for autosomes. 

For example, if you have a dataset split by chromosome in the following way: 

- `my_dataset_chr1.bgen`
- `my_dataset_chr2.bgen`
- `my_dataset_chr3.bgen`
- ...

You can specify the input dataset as `my_dataset_chr{CHROM}.bgen` and the pipeline will automatically replace `{CHROM}` with the chromosome number.

## Performance tips 

Note that the pipeline will do highly parallelized access to the input dataset. Thus, especially if you are using a single file as input and when step 2 chunk size is small, it is reccomended to use a high performance filesystem optimized for parallel access.

Regenie step2 will run faster on bgen v1.2 with 8 bits encoding. You can convert existing data using plink2 with `--export bgen-1.2 'bits=8'` option