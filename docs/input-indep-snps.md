# QCed independent SNPs

To perform regenie step 1 you should provide a {bed,bim,fam} dataset containing independent SNPs with reasonable QC. The step 1 input dataset is defined by `genotypes_array` parameter. This need to be specified as prefix, so if you a dataset named `my_dataset.bed/bim/fam` you should specify `my_dataset` as input.

Ideally, this would contain ~500k SNPs obtained after standard QC and pruning and MUST contain **less than 1M SNPs**. Usual pre-processing steps are:

- SNP missing rate < 0.05
- sample missing rate < 0.1
- MAF > 0.01
- HWE p-value > 1e-15
- variants pruning for LD (e.g. `--indep-pairwise 1000 100 0.8` in plink2)

An additional QC will be automatically performed on this file since step1 requires strict filtering criteria. Additionally, pruning can also be performed by setting `prune_enabled = true`.
