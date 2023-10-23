# nf-pipeline-regenie

A nextflow pipeline to perform genome-wide association studies (GWAS) and rare variant association analysis using [regenie](https://github.com/rgcgithub/regenie) at high speed.

<img src="images/regenie_pipeline.png" width="800" height="300">

## Main features

- The pipeline allows great flexibility in how operation are parallelized to allow balacing resource usage and speed and to accomodate different dataset sizes. The default settings are optimized for biobank scale data, but you can adjust resources as needed to further cut down run time or to process smaller datasets (see the [optimize parallelization section](parallelization.md)).
- All major data types are accepted as input, including plink1 binary dataset (bed/bim/fam), plink2 binary dataset (pgen/pvar/psam), bgen format (bgen/bgi/sample), and vcf.gz format.
- The pipeline can perform both standard GWAS analysis on single variants, and aggregated rare variant tests using burden test and any of the tests available in regenie, namely skat, skato, sakto-acat, acatv, acato, acato-full.
- Results include summary statistics per phenotype, and also filtered tophits / loci annotated with nearby genes and an HTML report for each phenotype with Manhattan plot and regional plots for the best loci.
- You can run a [single analysis project](running-modes.md#single-project) that can combine any number of phenotypes and covariates in a single analysis as long as you use a single model and phenotypes are all quantitative or all binary.
- Alternatively you can configure multiple analysis at once using [**multi models mode**](running-modes.md#multi-models) or [**multi projects mode**](running-modes.md#multi-projects). Using the multi models / projects mode it is possible to fully automate the test of multiple association models for a cohort. In multi-model mode the pipeline will take care of setting up uniform analysis groups given the provided models, while multi-project mode allows to run multiple independent analysis in parallel given a general configuration table.

## How to use

The suggested way to run the pipeline is to create a config file defining your computations environment (see the [hpc profile section](hpc-profile.md)) and a config file for your project (see the [main parameters section](main-parameters.md)). You can use the templates provided in the `templates` folder.

Then you can invoke the pipeline using `nextflow run HTGenomeAnalysisUnit/nf-pipeline-regenie -profile singularity,myprofile -c your_project.conf -c your_profile.conf`

See the [quick start section](quick-start.md) for a minimal example.

## Credits

The original concept is based on this amazing [github repository](https://github.com/genepi/nf-gwas) from Institute of Genetic Epidemiology, Innsbruck maintained by Sebastian Schönherr and Lukas Forer.
