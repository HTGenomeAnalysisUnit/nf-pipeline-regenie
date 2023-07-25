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