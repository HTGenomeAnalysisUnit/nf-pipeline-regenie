# REGENIE GWAS Pipeline - WDL Conversion

WDL (Workflow Description Language) conversion of the nf-pipeline-regenie Nextflow pipeline (v1.9.4).

## Files

| File | Description |
|------|-------------|
| `regenie_gwas.wdl` | **Main workflow** - entry point, orchestrates the full pipeline |
| `structs.wdl` | Struct definitions (PhenoMeta, CovarMeta, GenotypeFiles, etc.) |
| `tasks.wdl` | All task definitions (~25 tasks mapping to Nextflow processes) |
| `regenie_step1.wdl` | Subworkflow: QC, pruning, regenie step 1 (split L0/L1) |
| `regenie_step2_gwas.wdl` | Subworkflow: GWAS step 2 with scatter over chromosomes/chunks |
| `regenie_step2_rarevars.wdl` | Subworkflow: Rare variant step 2 with scatter over gene chunks |
| `process_gwas_results.wdl` | Subworkflow: Filter, annotate, and clump GWAS results |
| `process_rarevar_results.wdl` | Subworkflow: Filter rare variant results + corrected p-values |
| `inputs.json` | Template inputs JSON with all parameters |
| `options.json` | Runtime options for Cromwell/other WDL engines |
| `backends/` | Backend configuration files for different executors |
| `backends/local.conf` | Cromwell config: Local execution (Docker) |
| `backends/slurm.conf` | Cromwell config: SLURM cluster (Singularity) |
| `backends/lsf.conf` | Cromwell config: LSF cluster (Singularity) |
| `backends/miniwdl.cfg` | miniwdl config: SLURM via miniwdl-slurm plugin |

## Architecture Mapping

### Nextflow → WDL Correspondence

| Nextflow Component | WDL Equivalent |
|---|---|
| `main.nf` | `regenie_gwas.wdl` (main workflow) |
| `workflows/prepare_project.nf` | Input validation section in `regenie_gwas.wdl` |
| `workflows/variant_analysis.nf` | Step 2 sections in `regenie_gwas.wdl` |
| `subworkflow/regenie_step1.nf` | `regenie_step1.wdl` |
| `subworkflow/regenie_step2.nf` | `regenie_step2_gwas.wdl` + `regenie_step2_rarevars.wdl` |
| `subworkflow/prepare_step2_data.nf` | Handled via `GenotypeFiles` struct inputs |
| `subworkflow/split_data.nf` | `MakeVariantsChunks` / `MakeGenesChunks` tasks |
| `subworkflow/process_results.nf` | `process_gwas_results.wdl` + `process_rarevar_results.wdl` |
| `subworkflow/clump_results.nf` | Clumping tasks inside `process_gwas_results.wdl` |
| `modules/local/*.nf` | Individual tasks in `tasks.wdl` |

### Key Design Differences

1. **Channel operations → Scatter/Gather**: Nextflow channels with `.map`, `.combine`, `.groupTuple`, `.branch` are replaced by WDL `scatter` blocks with explicit array indexing.

2. **Dynamic fan-out**: Nextflow dynamically creates channels from file patterns. In WDL, genotype files must be provided as structured `GenotypeFiles` arrays in the inputs JSON.

3. **Multi-model/Multi-project modes**: The Nextflow pipeline has 3 execution modes (single, models_table, projects_table). The WDL version supports single-project execution. For multi-project, invoke the workflow multiple times.

4. **Conditional execution**: Nextflow `if (params.X)` blocks map to WDL `if (defined(X))` conditional blocks.

5. **Docker containers**: The pipeline uses 3 containers:
   - `htgenomeanalysisunit/gwas-nf-pipeline:0.7` (plink2, bedtools, python scripts)
   - `ghcr.io/rgcgithub/regenie/regenie:v4.0.gz` (regenie)
   - `edg1983/gwas-rmd-reports:v3.0` (quarto reports)

6. **bin/ scripts**: Python/R scripts from `bin/` must be either:
   - Included in the Docker container (recommended)
   - Passed as `File` inputs to the workflow

## Running

### With Cromwell — Local Backend
```bash
java -Dconfig.file=backends/local.conf \
    -jar cromwell.jar run regenie_gwas.wdl \
    -i inputs.json \
    -o options.json
```

### With Cromwell — SLURM Backend
```bash
# Tasks are submitted to SLURM via sbatch.
# Docker images are auto-converted to Singularity.
java -Dconfig.file=backends/slurm.conf \
    -jar cromwell.jar run regenie_gwas.wdl \
    -i inputs.json \
    -o options.json
```

To customize SLURM settings, either edit `backends/slurm.conf` or use runtime
attributes in the inputs JSON:
```json
{
    "RegenieGWAS.RegenieStep1.runtime_attr_override": {
        "slurm_partition": "gpu",
        "slurm_account": "myproject",
        "time_minutes": 1440
    }
}
```

### With Cromwell — LSF Backend
```bash
# Tasks are submitted to LSF via bsub.
# Docker images are auto-converted to Singularity.
java -Dconfig.file=backends/lsf.conf \
    -jar cromwell.jar run regenie_gwas.wdl \
    -i inputs.json \
    -o options.json
```

### With miniwdl — Local
```bash
miniwdl run regenie_gwas.wdl \
    -i inputs.json
```

### With miniwdl — SLURM
```bash
# Requires: pip install miniwdl-slurm
# Uses backends/miniwdl.cfg for configuration
MINIWDL__CFG=backends/miniwdl.cfg miniwdl run regenie_gwas.wdl \
    -i inputs.json
```

### On DNAnexus (dxCompiler)
```bash
java -jar dxCompiler.jar compile regenie_gwas.wdl \
    -project project-xxx \
    -folder /pipelines/ \
    -inputs inputs.json
```

## Executor Backends

### Overview

| Backend | Engine   | Container Tech | Config File             | Cluster Required |
|---------|----------|----------------|-------------------------|------------------|
| Local   | Cromwell | Docker         | `backends/local.conf`   | No               |
| SLURM   | Cromwell | Singularity    | `backends/slurm.conf`   | Yes              |
| LSF     | Cromwell | Singularity    | `backends/lsf.conf`     | Yes              |
| Local   | miniwdl  | Docker         | (built-in)              | No               |
| SLURM   | miniwdl  | Singularity    | `backends/miniwdl.cfg`  | Yes              |

### SLURM Configuration (Cromwell)

The SLURM backend submits tasks via `sbatch` and monitors them with `squeue`.
Docker images are automatically pulled and converted to Singularity SIF files.

Key settings in `backends/slurm.conf`:
- **`concurrent-job-limit`**: Max simultaneous SLURM jobs (default: 50)
- **`runtime-attributes`**: Optional per-task settings:
  - `slurm_partition` — target partition
  - `slurm_account` — billing account
  - `slurm_qos` — quality-of-service class
  - `slurm_extra` — any additional sbatch flags
  - `time_minutes` — wall time limit (default: 480 = 8h)

### LSF Configuration (Cromwell)

The LSF backend submits tasks via `bsub` and monitors them with `bjobs`.

Key settings in `backends/lsf.conf`:
- **`concurrent-job-limit`**: Max simultaneous LSF jobs (default: 50)
- **`runtime-attributes`**: Optional per-task settings:
  - `lsf_queue` — target queue
  - `lsf_project` — billing project
  - `lsf_extra` — any additional bsub flags
  - `time_minutes` — wall time limit (default: 480 = 8h)

### SLURM Configuration (miniwdl)

Requires the `miniwdl-slurm` plugin:
```bash
pip install miniwdl miniwdl-slurm
```

Edit `backends/miniwdl.cfg` to set:
- **`[singularity] image_cache`**: Directory for cached Singularity images
- **`[singularity] extra_args`**: Bind mounts for your cluster (e.g., `--bind /scratch:/scratch`)
- **`[slurm] extra_args`**: Default sbatch flags (e.g., `--partition=normal`)

### Singularity Image Caching

On HPC clusters (SLURM/LSF), Docker images are automatically converted to
Singularity SIF files. These are cached to avoid repeated pulls:

- **Cromwell**: Images cached in `${cwd}/singularity_images/` per task execution
- **miniwdl**: Images cached in the directory set by `[singularity] image_cache`

For shared cluster environments, consider setting a shared cache directory:
```bash
# Cromwell: edit backends/slurm.conf or backends/lsf.conf
SINGULARITY_IMAGE_DIR="/shared/singularity_cache"

# miniwdl: edit backends/miniwdl.cfg
[singularity]
image_cache = /shared/singularity_cache
```

## Important Notes

1. **This is a best-effort conversion** - the WDL has not been validated against the original Nextflow test suite. Testing is required before production use.

2. **Data preparation**: The Nextflow pipeline auto-discovers files via glob patterns (e.g., `*.{bed,bim,fam}`). In WDL, you must explicitly list all input files in the `GenotypeFiles` struct array.

3. **Results concatenation**: The Nextflow pipeline groups results by phenotype using `groupTuple`. The WDL version passes all results to `ConcatStep2Results` which filters by phenotype name in filenames. This may need refinement for complex multi-chromosome setups.

4. **Local executor processes**: Nextflow `CHECK_PROJECT` and `CHECK_CHANNEL_SIZE` run with `executor 'local'` for validation. These are omitted in WDL as they have no command - their validation logic should be implemented in the calling code or as pre-flight checks.
