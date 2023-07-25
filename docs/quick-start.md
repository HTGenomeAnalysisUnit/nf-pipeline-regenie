# Quick Start

1. Create a folder for your project (e.g. `yourproject`) 

2. Prepare and configure the required [input data for step 2](input-full-data.md), usually an imputed or sequencing dataset, and [step 1](input-indep-snps.md), usually a QCed and pruned dataset.

3. Prepare the [config files](main-parameters.md) for your project in your project folder, including a config file to [define the profile](hpc-profile.md) for your computational environment. You can use the templates provided in the `templates` folder.

4. Invoke the pipeline

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