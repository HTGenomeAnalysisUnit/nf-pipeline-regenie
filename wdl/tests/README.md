# Running WDL Tests

This directory contains test input JSON files and a runner script that mirror the
Nextflow `run-tests.sh` test suite.

## Test Coverage

### Single-project mode (Mode 1)

| Test   | Description                                          | Scenarios                                    |
|--------|------------------------------------------------------|----------------------------------------------|
| test1  | BGEN input, quant phenotype, 3 covars                | Split step2, clumping, report                |
| test2  | BGEN input, binary phenotype, 3 covars               | Split step2, clumping, report                |
| test3  | BGEN input, quant phenotype, no covariates           | Split step2, clumping                        |
| test4  | BGEN input **by chromosome**, quant phenotype        | Split step2, per-chr genotype files (\*)     |
| test5  | BGEN input, quant phenotype                          | **No split** step2                           |
| test6  | BGEN input **by chromosome**, quant phenotype        | **No split** step2, per-chr genotype (\*)    |
| test7  | **PGEN** input, quant phenotype                      | Split step2                                  |
| test8  | **VCF** input, quant phenotype                       | Split step2                                  |
| test9  | **VCF** input **by chromosome**, quant phenotype     | Split step2, per-chr VCF files (\*)          |
| test10 | BGEN input, quant phenotype                          | **Premade step1 predictions**                |
| test11 | BGEN input, quant phenotype                          | **Skip step1 predictions entirely**          |
| test12 | PGEN input **by chromosome**, quant phenotype        | Split step2, per-chr genotype files          |
| test13 | PGEN input **multiple chunk files**, quant phenotype | Split step2, 5 chunk files                   |
| test14 | BGEN input **multiple chunk files**, quant phenotype | Split step2, 5 chunk files (\*)              |
| test15 | BGEN input **multiple chunk files**, quant phenotype | **No split** step2, 5 chunk files (\*)       |
| test16 | BGEN input, quant phenotype, 3 covars                | Split step2, **no clumping**, report         |

(\*) Requires index files that may not ship with the test data. See [Prerequisites](#bgen--vcf-index-files) below.

### Models table mode (Mode 2)

| Test   | Description                                          | Scenarios                                    |
|--------|------------------------------------------------------|----------------------------------------------|
| test20 | Models table with 5 models (quant + binary)          | Split step2, uses `models_table`             |

### Projects array mode (Mode 3)

| Test   | Description                                          | Scenarios                                    |
|--------|------------------------------------------------------|----------------------------------------------|
| test21 | Two projects (quant + binary) via `ProjectConfig`    | Split step2, uses `projects` array           |

### Tests not included (require additional handling)

| Test   | Reason                                                              |
|--------|---------------------------------------------------------------------|
| test17 | BGEN no BGI — WDL struct requires .bgi; pre-index with `bgenix`    |
| test18 | BGEN no sample file + double-ID — WDL struct requires sample file  |
| test19 | BGEN no sample + no BGI — combines test17 + test18 issues          |

## Prerequisites

### Option A: Cromwell
```bash
# Download Cromwell JAR
wget https://github.com/broadinstitute/cromwell/releases/download/87/cromwell-87.jar
export CROMWELL_JAR=/path/to/cromwell-87.jar
```

### Option B: miniwdl
```bash
pip install miniwdl
```

### BGEN / VCF index files
Some tests (marked with \* above) reference per-chromosome or per-chunk BGEN `.bgi`
or VCF `.tbi` index files that the Nextflow pipeline generates on the fly but the
WDL pipeline requires upfront. Generate them before running those tests:
```bash
# BGEN index files (for tests 4, 6, 14, 15)
for f in tests/input/pipeline/example_chr*.bgen tests/input/pipeline/example_chunk*.bgen; do
    bgenix -g "$f" -index
done

# VCF index files (for test 9)
for f in tests/input/pipeline/example_chr*.vcf.gz; do
    tabix -p vcf "$f"
done
```

### Docker Images
The pipeline uses these Docker containers (must be available):
- `htgenomeanalysisunit/gwas-nf-pipeline:0.7` — plink2, bedtools, python scripts
- `ghcr.io/rgcgithub/regenie/regenie:v4.0.gz` — regenie
- `edg1983/gwas-rmd-reports:v3.0` — quarto reports

Pull them beforehand:
```bash
docker pull htgenomeanalysisunit/gwas-nf-pipeline:0.7
docker pull ghcr.io/rgcgithub/regenie/regenie:v4.0.gz
docker pull edg1983/gwas-rmd-reports:v3.0
```

## Running Tests

### Run all tests with the test runner
```bash
cd wdl/

# With Cromwell
./run-tests-wdl.sh cromwell

# With miniwdl
./run-tests-wdl.sh miniwdl
```

### Run a specific test
```bash
# Run only test1 with Cromwell
./run-tests-wdl.sh cromwell local test1

# Run only test7 with miniwdl
./run-tests-wdl.sh miniwdl local test7

# Run models-table test (Mode 2)
./run-tests-wdl.sh cromwell local test20

# Run projects-array test (Mode 3)
./run-tests-wdl.sh cromwell local test21
```

### Run a single test manually

**With Cromwell:**
```bash
cd wdl/

# The runner resolves relative paths automatically. For manual runs:
java -jar cromwell.jar run regenie_gwas.wdl \
  --inputs tests/test1_inputs.json \
  --options options.json
```

**With miniwdl:**
```bash
cd wdl/

miniwdl run regenie_gwas.wdl \
  --input tests/test1_inputs.json \
  --dir output/test1
```

## Input JSON Structure

Each test JSON maps Nextflow parameters to WDL inputs. Key differences from Nextflow:

| Nextflow Parameter                    | WDL Input                           | Notes                                      |
|---------------------------------------|-------------------------------------|--------------------------------------------|
| `genotypes_array = "prefix"`          | `genotypes_array_bed/bim/fam`       | Split into 3 explicit file paths           |
| `genotypes_imputed = "path"`          | `genotypes_imputed` (array)         | Array of `GenotypeFiles` structs           |
| `genotypes_imputed = "path_{CHROM}"`  | Multiple entries in array           | One `GenotypeFiles` struct per chromosome  |
| `chromosomes = "1,2"`                | `chromosomes = ["1","2"]`           | String → Array[String]                     |
| `covariates_filename = 'NO_COV_FILE'`| Omit `covariates_filename`          | Simply don't include the optional input    |
| `regenie_premade_predictions = glob`  | `premade_prediction_files` array    | Explicit file list instead of glob pattern |

### GenotypeFiles Struct
```json
{
    "file_prefix": "example",
    "primary_file": "/path/to/example.bgen",
    "secondary_file": "/path/to/example.bgen.bgi",
    "tertiary_file": "/path/to/example.sample",
    "chromosome": "ONE_FILE"
}
```

For **BGEN**: primary=.bgen, secondary=.bgen.bgi, tertiary=.sample  
For **PGEN**: primary=.pgen, secondary=.pvar, tertiary=.psam  
For **VCF**:  primary=.vcf.gz, secondary=.vcf.gz.tbi, tertiary=.sample  
For **BED**:  primary=.bed, secondary=.bim, tertiary=.fam  

Set `chromosome` to `"ONE_FILE"` for single whole-genome files, or `"1"`, `"2"`, etc. for per-chromosome files.

## Path Resolution

The test input JSONs use **relative paths** (e.g., `../tests/input/pipeline/...`).
The `run-tests-wdl.sh` script automatically resolves these to absolute paths before
passing them to the WDL engine. If running manually, ensure paths are absolute.
