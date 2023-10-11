# LD panel (optional)

If you are analysing a large dataset with more than 50k samples, to speed up LD computation for clumping we suggest to prepare by chromosome bed/bim/fam files with the same variants present in the full genotype data but only a subset of samples. Then you can specify a pattern to these files using the `ld_panel` parameter like `/ld_panel/chr{CHROM}_ld`.

The `{CHROM}` is automatically substituted with number 1-23 when the pipeline is running. This is only processed when the `clumping` option is active. If the `ld_panel` parameter is not set and clumping is active the pipeline will use the full genotype data to estimate LD. Note that this will result in very long run time for huge datasets so providing an LD panel is highly reccomended when sample size is above 50k.

You can generate LD files using plink2 and the `--keep` option to extract a subset of unrelated samples from the full dataset. For example:

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

## When LD panel is not available

If you do not provide an LD panel, the pipeline will use the full input dataset to compute LD. First, input file(s) will be converted to bed format using plink2. Then if the input files are organized by chromosomes using the `{CHROM}` placeholder, only the relevant files will be used to perform clumping. Otherwise, if there are multiple input files, all converted dataset will be merged into a single dataset and this is used to perform clumping.

In any case, clumping procedure is then performed parallelized by chromosome and resulting annotated loci are then merged into a single file per phenotype.
