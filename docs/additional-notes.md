
# Additional notes

## MAC filtering

Note that the pipeline impose a MAC (100) and a MAF (0.01) filter on genotyped data used for step1, as reccomended by [regenie UKBB tutorial](https://rgcgithub.github.io/regenie/recommendations/#step-1), and a min MAC 50 at step2 as suggested in the recent large proteome analysis of UKBB (see [the preprint](https://www.biorxiv.org/content/10.1101/2022.06.17.496443v1)). 

The min MAC for variants accepted in step2 can be adjusted by setting the `regenie_rarevar_min_mac` and `regenie_gwas_min_mac` parameters in your profile.

This MAC filters are applied after subsetting the data to contain only samples seen in the phenotype table. FID/IID are extracted from the pheno file and used with `--keep` in genotyped data QC to generate a QCed dataset used in regenie step 1. Then in step 2, regenie automatically removes any sample not seen in step 1 predictions before applying the MAC filter.

As a result of this process, if you use only a subset of samples present in the original imputed data, it is possible that some variants will be filtered out due to low MAC. Check log files and the results table to see exactly how many variants were considered in the analysis.

## Clumping

When clumping is active, the pipeline will save clumped data and resulting loci with genes annotation in the `toploci` folder. Note that if there are multiple identical SNP IDs clumping will fail. To avoid this you can for example modify your SNP ids to include ref/alt alleles as follows: `[SNPID]_[A0]_[A1]`.

If you are using LD panel files as described in the [LD panel section](input-ld-panel.md), please ensure that SNP ids are concordant between bed/bim/fam files in the LD_panel and the input  dataset for GWAS analysis, otherwise SNPs that are not found will be dropped from clumping report.