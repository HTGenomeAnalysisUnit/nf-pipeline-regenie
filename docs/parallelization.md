# Parallelization

The pipeline optimize the run time by splitting heavy tasks into small chunks. The default settings are optimized for biobank scale data, but you can adjust resources as needed to further cut down run time or to process smaller datasets. You can change the default minimum resources used at each step by editing the file `conf/base.config`, or providing a new config file with the `-c` option when you run the pipeline.

First of all, keep in mind that the analysis of binary phenotypes rquires much more resources and take longer times. In general we suggest to analyze not more than 100 quantitative phenotypes or 50 binary phenotypes in a single run. If you need to analyze more phenotypes, you can take advatnge of the multi-model run mode to split the analysis in multiple runs.

You can optimize the pipeline scaling by adjusting the following parameters:

- `step1_n_chunks`: this controls the number of chunks for step1 L0 regression. The default of 100 is usually fine in most case and can eventually be reduced for small datasets with less than 10 thousands samples.
- step1 L1 is then performed by phenotype. This behaviour can not be changed. If you are analyzing a large number of binary phenotypes (more than 50 in a single run) on biobank scale date (more than 200k individuals) you will probably need to increase the min memory for `step1_runl1` label to 16 or 24G.
- `step2_gwas_chunk_size`: this controls the number of variants analyzed in each chunk during regenie step2. Reducing this number increase the number of concurrent chunks and reduce run time per chunk, thus increasing overall performance when enough resources are available. Default is 50000. This can be increased to 100000 or even 250000 for small datasets with less than 100k individuals.
- `step2_rarevar_chunk_size`: this controls the number of genes from the input snpset file analyzed in each chunk during regenie step2 for rare variant analysis. Reducing this number increase the number of concurrent chunks and reduce run time when enough resources are available. Default is 100.

Another important aspect to consider is the maximum number of jobs submitted simultaneously. This can be asjuted using the `queueSize` parameter in the executor scope when you create your custom profile. If you have a large cluster with many resources available, you can set this to a large value to submit all jobs simultaneously. This will reduce the overall run time when enough resources are available.

To further increase speed, you can eventually adjust the default amount of resources requested for regenie step1 and regenie step2 operations by editing the file `conf/base.config`. The default here are sensible for most cases, but you can increase the amount of memory and CPUs requested for each task to reduce computation time. In particular you can act on the following labels in the configuration file to adjust the resources used as starting point in step1 and step2: step1_runl0, step1_runl1, step2_gwas, step2_rarevar. Resources allocation is dynamic, so when a process fail due to OOM or OOT it will be retried automatically increasing the resources.

In our tests using the default settings one can analyze 100 quantitative phenotypes on the full UKBB dataset (about 500k individuals and 45M SNPs) including loci identification, but withouth HTML reports, in ~7h under normal HPC load considering a limit of 200 concurrent tasks (peak of 800 CPUs usage).
