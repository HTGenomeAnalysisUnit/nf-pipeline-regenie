# Parallelization

The pipeline optimize the run time by splitting heavy tasks into small chunks. You can control the exact beahviour of this process by adjusting the following parameters:

- `step1_n_chunks`: this controls the number of chunks for step1 L0 regression. The default of 100 is usually fine in most case and can eventually be reduced for small datasets with few thousands samples.
- `step2_gwas_chunk_size`: this controls the number of variants analyzed in each chunk during regenie step2. Reducing this number increase the number of concurrent chunks and reduce run time when enough resources are available. Default is 1000000.
- `step2_rarevar_chunk_size`: this controls the number of genes from the input snpset file analyzed in each chunk during regenie step2 for rare variant analysis. Reducing this number increase the number of concurrent chunks and reduce run time when enough resources are available. Default is 200.

You can also control the number of jobs submitted simultaneously by adjusting the `queueSize` parameter in the executor scope when you create your custom profile. If you have a large cluster with many resources available, you can set this toa large value to submit all jobs simultaneously. This will reduce the overall run time when enough resources are available.

To further increase speed, you can eventually adjust the default amount of resources requested for regenie step1 and regenie step2 operations by editing the file `conf/base.config`. The default here are sensible for most cases, but you can increase the amount of memory and CPUs requested for each task to reduce computation time.

In our tests using the default settings one can analyze 100 quantitative phenotypes on the full UKBB dataset (about 500k individuals and 40M SNPs) in ~3h under normal HPC load considering a limit of 200 concurrent tasks. 