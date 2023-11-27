# Quick Start

1. Create a folder for your project (e.g. `yourproject`)

2. Prepare a tab-separated table of phenotypes and eventually covariates (see the [input section](input-phenotype-file.md)).

3. Prepare and configure the required [input data for step 2](input-full-data.md), usually an imputed or sequencing dataset, and [step 1](input-indep-snps.md), usually a QCed and pruned dataset. You can eventually prepare also a [set of files for LD computation](input-ld-panel.md), suggested if you want to run loci clumping when analyzing a large dataset with > 100k samples.

4. If you want to perform a multi-models or multi-projects execution, prepare the [models table](input-models-table.md) or [projects table](input-projects-table.md) to describe your analyses.

5. If you want to perform conditional or interaction analysis, prepare the additional genotype dataset and eventually conditional variant list as described in the [conditional analysis section](input-conditional-analysis.md).

6. Prepare the necessary config files, using the templates provided in the `templates` folder:
   1. A [config file](main-parameters.md) describing settings and inputs for your project.
   2. A config file to [define the profile](hpc-profile.md) for your computational environment.
   3. Optionally, you can also add configuration to enable execution monitoring using [Nextflow Tower](tower-monitoring.md)

7. Invoke the pipeline using `nextflow run HTGenomeAnalysisUnit/nf-pipeline-regenie`

Usually, you want to prepare a script to submit the pipeline in your project folder. In this example we use `sbatch` submission system, but this can be adapted to any scheduler. `myprofile` corresponds to a profile you created for your computational environment:

```bash
#!/bin/bash
#SBATCH --job-name nf-regenie
#SBATCH --output nf-regenie_master_%A.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 1-00:00:00

module load nextflow/22.10.1 singularity/3.8.5

export NXF_OPTS="-Xms1G -Xmx8G" 
nextflow run HTGenomeAnalysisUnit/nf-pipeline-regenie \
   -profile singularity,myprofile -c your_project.conf
```

Alternatively, you can clone the latest pipeline version using

`git clone --depth 1 https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie.git`

This will create a new folder called `nf-pipeline-regenie` in the current folder containing all the pipeline files.

You can eventually chose a specific version of the pipeline using the `--branch` option

`git clone --depth 1 --branch v1.8 https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie.git`