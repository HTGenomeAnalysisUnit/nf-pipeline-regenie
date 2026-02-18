version 1.0

import "structs.wdl" as Structs

## ============================================================================
## Task Definitions for the REGENIE GWAS Pipeline
## ============================================================================

## ---------- Input Validation Tasks ----------

task ValidatePhenotypes {
    input {
        String project_id
        File phenotypes_file
        PhenoMeta pheno_meta

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 20,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    String out_prefix = basename(phenotypes_file, ".txt")

    command <<<
        set -euo pipefail
        RegenieValidateInput.py \
            --input ~{phenotypes_file} \
            --output ~{out_prefix}.pheno.validated.txt \
            --type phenotype
    >>>

    output {
        File validated_phenotypes = "~{out_prefix}.pheno.validated.txt"
        File validation_log = "~{out_prefix}.pheno.validated.log"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task ValidateCovariates {
    input {
        String project_id
        File covariates_file
        CovarMeta covar_meta
        AccessoryFiles accessory_files

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 20,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    String out_prefix = basename(covariates_file, ".txt")

    command <<<
        set -euo pipefail
        RegenieValidateInput.py \
            --input ~{covariates_file} \
            --output ~{out_prefix}.cov.validated.txt \
            --type covariate
    >>>

    output {
        File validated_covariates = "~{out_prefix}.cov.validated.txt"
        File validation_log = "~{out_prefix}.cov.validated.log"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task SetupMultipleRuns {
    input {
        File pheno_chunker_r
        File prepare_projects_py
        File traits_table
        File models_table
        File fam_file
        Int pheno_chunk_size = 50
        Float missing_tolerance = 0.1

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 20,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    command <<<
        set -euo pipefail

        # Step 1: Run phenotype chunker (R) and project builder (Python)
        Rscript --vanilla ~{pheno_chunker_r} \
            --file ~{traits_table} \
            --model_file ~{models_table} \
            --sample_file ~{fam_file} \
            -s ~{pheno_chunk_size} \
            -t ~{missing_tolerance}

        python ~{prepare_projects_py} master_table.tsv

        # Step 2: Post-process for WDL consumption
        # Reorganize outputs into parallel arrays with deterministic ordering
        python3 <<'PYEOF'
import csv
import os
import shutil

projects = []
with open('analysis.conf') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        projects.append(row)

n = len(projects)

# Write parallel metadata arrays (one value per line)
with open('project_ids.txt', 'w') as f:
    for p in projects:
        f.write(p['project_id'] + '\n')

with open('pheno_cols.txt', 'w') as f:
    for p in projects:
        f.write(p['pheno_cols'] + '\n')

with open('pheno_binary.txt', 'w') as f:
    for p in projects:
        f.write(p['pheno_binary'] + '\n')

with open('pheno_model.txt', 'w') as f:
    for p in projects:
        f.write(p['pheno_model'] + '\n')

with open('cov_cols.txt', 'w') as f:
    for p in projects:
        f.write(p.get('cov_cols', '') + '\n')

with open('cov_cat_cols.txt', 'w') as f:
    for p in projects:
        f.write(p.get('cov_cat_cols', '') + '\n')

# Copy phenotype and covariate files with deterministic zero-padded names
# This ensures glob() returns them in the same order as the metadata arrays
has_cov_lines = []
for i, p in enumerate(projects):
    # Copy phenotype file
    shutil.copy(p['pheno_file'], f'phenodata_{i:04d}.tsv')

    # Handle covariate file
    cov_file = p.get('cov_file', 'NO_COV_FILE')
    if cov_file not in ('NO_COV_FILE', 'NA', ''):
        shutil.copy(cov_file, f'covardata_{i:04d}.tsv')
        has_cov_lines.append('true')
    else:
        # Create empty marker file so glob count matches
        with open(f'covardata_{i:04d}.tsv', 'w') as cf:
            cf.write('')
        has_cov_lines.append('false')

with open('has_cov.txt', 'w') as f:
    for line in has_cov_lines:
        f.write(line + '\n')

print(f'Prepared {n} project configurations for WDL scatter')
PYEOF
    >>>

    output {
        File analysis_config = "analysis.conf"
        File master_table = "master_table.tsv"
        Array[String] project_ids = read_lines("project_ids.txt")
        Array[String] pheno_cols_list = read_lines("pheno_cols.txt")
        Array[String] pheno_binary_list = read_lines("pheno_binary.txt")
        Array[String] pheno_model_list = read_lines("pheno_model.txt")
        Array[File] pheno_files = glob("phenodata_*.tsv")
        Array[String] cov_cols_list = read_lines("cov_cols.txt")
        Array[String] cov_cat_cols_list = read_lines("cov_cat_cols.txt")
        Array[File] cov_files = glob("covardata_*.tsv")
        Array[String] has_cov_list = read_lines("has_cov.txt")
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

## ---------- Step 1 Pre-processing Tasks ----------

task QcFilterGenotyped {
    input {
        String project_id
        File bed_file
        File bim_file
        File fam_file
        File phenos_tsv

        String qc_maf = "0.01"
        String qc_mac = "100"
        String qc_geno = "0.05"
        String qc_hwe = "1e-15"
        String qc_mind = "0.1"

        RuntimeAttributes runtime_attr = object {
            cpu: 4,
            memory: "8 GB",
            disk_size_gb: 100,
            max_retries: 3,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    String bed_basename = basename(bed_file, ".bed")

    command <<<
        set -euo pipefail

        # Symlink plink files to have consistent naming
        ln -sf ~{bed_file} ~{bed_basename}.bed
        ln -sf ~{bim_file} ~{bed_basename}.bim
        ln -sf ~{fam_file} ~{bed_basename}.fam

        grep -v "FID" ~{phenos_tsv} | cut -f1,2 > samples.list

        plink2 \
            --bfile ~{bed_basename} \
            --keep samples.list \
            --maf ~{qc_maf} \
            --mac ~{qc_mac} \
            --geno ~{qc_geno} \
            --hwe ~{qc_hwe} \
            --mind ~{qc_mind} \
            --write-snplist --write-samples --no-id-header \
            --out ~{bed_basename}.qc \
            --make-bed \
            --threads ~{runtime_attr.cpu} \
            --memory ~{4000}
    >>>

    output {
        File qc_bed = "~{bed_basename}.qc.bed"
        File qc_bim = "~{bed_basename}.qc.bim"
        File qc_fam = "~{bed_basename}.qc.fam"
        File qc_log = "~{bed_basename}.qc.log"
        File qc_snplist = "~{bed_basename}.qc.snplist"
        File qc_id = "~{bed_basename}.qc.id"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task PruneGenotyped {
    input {
        String project_id
        File bed_file
        File bim_file
        File fam_file

        Float prune_maf = 0.01
        Int prune_window_kbsize = 1000
        Int prune_step_size = 100
        Float prune_r2_threshold = 0.8

        RuntimeAttributes runtime_attr = object {
            cpu: 4,
            memory: "8 GB",
            disk_size_gb: 100,
            max_retries: 3,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    String bim_basename = basename(bim_file, ".bim")

    command <<<
        set -euo pipefail

        ln -sf ~{bed_file} ~{bim_basename}.bed
        ln -sf ~{bim_file} ~{bim_basename}.bim
        ln -sf ~{fam_file} ~{bim_basename}.fam

        plink2 \
            --bfile ~{bim_basename} \
            --double-id --maf ~{prune_maf} \
            --indep-pairwise ~{prune_window_kbsize} ~{prune_step_size} ~{prune_r2_threshold} \
            --out ~{bim_basename} \
            --threads ~{runtime_attr.cpu} \
            --memory ~{4000}

        plink2 \
            --bfile ~{bim_basename} \
            --extract ~{bim_basename}.prune.in \
            --double-id \
            --make-bed \
            --out ~{bim_basename}.pruned \
            --threads ~{runtime_attr.cpu} \
            --memory ~{4000}
    >>>

    output {
        File pruned_bed = "~{bim_basename}.pruned.bed"
        File pruned_bim = "~{bim_basename}.pruned.bim"
        File pruned_fam = "~{bim_basename}.pruned.fam"
        File pruned_log = "~{bim_basename}.pruned.log"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

## ---------- Regenie Step 1 Tasks ----------

task RegenieSplitL0 {
    input {
        String project_id
        File phenotypes_file
        PhenoMeta pheno_meta
        File covariates_file
        CovarMeta covar_meta
        AccessoryFiles accessory_files
        File bed_file
        File bim_file
        File fam_file

        Int regenie_bsize_step1 = 1000
        Int step1_n_chunks = 100
        Boolean phenotypes_delete_missings = false
        Boolean regenie_ref_first = false
        Int maxCatLevels = 10
        String? additional_geno_format

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "4 GB",
            disk_size_gb: 100,
            max_retries: 3,
            docker: "ghcr.io/rgcgithub/regenie/regenie:v4.0.gz"
        }
    }

    String bed_basename = basename(bed_file, ".bed")

    command <<<
        set -euo pipefail

        ln -sf ~{bed_file} ~{bed_basename}.bed
        ln -sf ~{bim_file} ~{bed_basename}.bim
        ln -sf ~{fam_file} ~{bed_basename}.fam

        # Build optional arguments
        COVARIANTS=""
        if [ "~{basename(covariates_file)}" != "NO_COV_FILE" ]; then
            COVARIANTS="--covarFile ~{covariates_file} --covarColList ~{covar_meta.cols}"
        fi

        CAT_COVARIATES=""
        if [ -n "~{covar_meta.cat_cols}" ] && [ "~{covar_meta.cat_cols}" != "NA" ]; then
            CAT_COVARIATES="--catCovarList ~{covar_meta.cat_cols}"
        fi

        DELETE_MISSINGS=""
        if [ "~{phenotypes_delete_missings}" = "true" ]; then
            DELETE_MISSINGS="--strict"
        fi

        REF_FIRST=""
        if [ "~{regenie_ref_first}" = "true" ]; then
            REF_FIRST="--ref-first"
        fi

        CONDITION_LIST=""
        ~{if defined(accessory_files.condition_list) then 'CONDITION_LIST="--condition-list ' + select_first([accessory_files.condition_list, ""]) + '"' else ''}

        regenie \
            --step 1 \
            --bed ~{bed_basename} \
            --phenoFile ~{phenotypes_file} \
            --phenoColList ~{pheno_meta.cols} \
            $COVARIANTS \
            $CAT_COVARIATES \
            $DELETE_MISSINGS \
            $REF_FIRST \
            --maxCatLevels ~{maxCatLevels} \
            $CONDITION_LIST \
            --bsize ~{regenie_bsize_step1} \
            --split-l0 regenie_step1,~{step1_n_chunks} \
            --out regenie_step1_splitl0
    >>>

    output {
        File master_file = "regenie_step1.master"
        Array[File] snplists = glob("regenie_step1*.snplist")
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task RegenieRunL0 {
    input {
        Int job_number
        String project_id
        File phenotypes_file
        PhenoMeta pheno_meta
        File covariates_file
        CovarMeta covar_meta
        AccessoryFiles accessory_files
        File master_file
        Array[File] snplists
        File bed_file
        File bim_file
        File fam_file

        Int regenie_bsize_step1 = 1000
        Boolean phenotypes_delete_missings = false
        Boolean regenie_force_step1 = false
        Boolean regenie_ref_first = false
        Boolean use_loocv = false
        Int maxCatLevels = 10
        String? additional_geno_format

        RuntimeAttributes runtime_attr = object {
            cpu: 2,
            memory: "12 GB",
            disk_size_gb: 100,
            max_retries: 3,
            docker: "ghcr.io/rgcgithub/regenie/regenie:v4.0.gz"
        }
    }

    String bed_basename = basename(bed_file, ".bed")
    String master_prefix = basename(master_file, ".master")

    command <<<
        set -euo pipefail

        ln -sf ~{bed_file} ~{bed_basename}.bed
        ln -sf ~{bim_file} ~{bed_basename}.bim
        ln -sf ~{fam_file} ~{bed_basename}.fam

        # Stage snplists and master file in current dir
        for f in ~{sep=" " snplists}; do
            ln -sf "$f" .
        done
        ln -sf ~{master_file} .

        COVARIANTS=""
        if [ "~{basename(covariates_file)}" != "NO_COV_FILE" ]; then
            COVARIANTS="--covarFile ~{covariates_file} --covarColList ~{covar_meta.cols}"
        fi

        CAT_COVARIATES=""
        if [ -n "~{covar_meta.cat_cols}" ] && [ "~{covar_meta.cat_cols}" != "NA" ]; then
            CAT_COVARIATES="--catCovarList ~{covar_meta.cat_cols}"
        fi

        DELETE_MISSINGS=""
        if [ "~{phenotypes_delete_missings}" = "true" ]; then
            DELETE_MISSINGS="--strict"
        fi

        FORCE_STEP1=""
        if [ "~{regenie_force_step1}" = "true" ]; then
            FORCE_STEP1="--force-step1"
        fi

        REF_FIRST=""
        if [ "~{regenie_ref_first}" = "true" ]; then
            REF_FIRST="--ref-first"
        fi

        USE_LOOCV=""
        if [ "~{use_loocv}" = "true" ]; then
            USE_LOOCV="--loocv"
        fi

        BINARY=""
        if [ "~{pheno_meta.binary}" = "true" ]; then
            BINARY="--bt"
        fi

        regenie \
            --step 1 \
            --bed ~{bed_basename} \
            --phenoFile ~{phenotypes_file} \
            --phenoColList ~{pheno_meta.cols} \
            $COVARIANTS \
            $CAT_COVARIATES \
            $DELETE_MISSINGS \
            $FORCE_STEP1 \
            $REF_FIRST \
            $USE_LOOCV \
            --maxCatLevels ~{maxCatLevels} \
            $BINARY \
            --threads ~{runtime_attr.cpu} \
            --bsize ~{regenie_bsize_step1} \
            --run-l0 ~{basename(master_file)},~{job_number} \
            --out regenie_step1_run-l0_~{job_number}
    >>>

    output {
        Array[File] l0_output = glob("~{master_prefix}_job~{job_number}_l0_*")
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task RegenieRunL1 {
    input {
        String project_id
        String single_pheno
        File phenotypes_file
        PhenoMeta pheno_meta
        File covariates_file
        CovarMeta covar_meta
        AccessoryFiles accessory_files
        File master_file
        Array[File] snplists
        Array[File] runl0_files
        File bed_file
        File bim_file
        File fam_file

        Int regenie_bsize_step1 = 1000
        Int niter = 30
        Boolean phenotypes_delete_missings = false
        Boolean regenie_force_step1 = false
        Boolean regenie_ref_first = false
        Boolean use_loocv = false
        Int maxCatLevels = 10
        String? additional_geno_format

        RuntimeAttributes runtime_attr = object {
            cpu: 6,
            memory: "12 GB",
            disk_size_gb: 100,
            max_retries: 3,
            docker: "ghcr.io/rgcgithub/regenie/regenie:v4.0.gz"
        }
    }

    String bed_basename = basename(bed_file, ".bed")

    command <<<
        set -euo pipefail

        ln -sf ~{bed_file} ~{bed_basename}.bed
        ln -sf ~{bim_file} ~{bed_basename}.bim
        ln -sf ~{fam_file} ~{bed_basename}.fam

        # Stage all L0 output, snplists, and master file in current dir
        for f in ~{sep=" " runl0_files}; do
            ln -sf "$f" .
        done
        for f in ~{sep=" " snplists}; do
            ln -sf "$f" .
        done
        ln -sf ~{master_file} .

        COVARIANTS=""
        if [ "~{basename(covariates_file)}" != "NO_COV_FILE" ]; then
            COVARIANTS="--covarFile ~{covariates_file} --covarColList ~{covar_meta.cols}"
        fi

        CAT_COVARIATES=""
        if [ -n "~{covar_meta.cat_cols}" ] && [ "~{covar_meta.cat_cols}" != "NA" ]; then
            CAT_COVARIATES="--catCovarList ~{covar_meta.cat_cols}"
        fi

        DELETE_MISSINGS=""
        if [ "~{phenotypes_delete_missings}" = "true" ]; then
            DELETE_MISSINGS="--strict"
        fi

        FORCE_STEP1=""
        if [ "~{regenie_force_step1}" = "true" ]; then
            FORCE_STEP1="--force-step1"
        fi

        REF_FIRST=""
        if [ "~{regenie_ref_first}" = "true" ]; then
            REF_FIRST="--ref-first"
        fi

        USE_LOOCV=""
        if [ "~{use_loocv}" = "true" ]; then
            USE_LOOCV="--loocv"
        fi

        BINARY=""
        if [ "~{pheno_meta.binary}" = "true" ]; then
            BINARY="--bt"
        fi

        regenie \
            --step 1 \
            --bed ~{bed_basename} \
            --phenoFile ~{phenotypes_file} \
            --phenoColList ~{pheno_meta.cols} \
            --l1-phenoList ~{single_pheno} \
            $COVARIANTS \
            $CAT_COVARIATES \
            $DELETE_MISSINGS \
            $FORCE_STEP1 \
            $REF_FIRST \
            $USE_LOOCV \
            --maxCatLevels ~{maxCatLevels} \
            $BINARY \
            --threads ~{runtime_attr.cpu} \
            --bsize ~{regenie_bsize_step1} \
            --niter ~{niter} \
            --run-l1 ~{basename(master_file)} \
            --keep-l0 --gz --verbose \
            --out regenie_step1_out_~{single_pheno}

        # Fix pred.list paths to be relative
        sed -i "s|$PWD/||g" regenie_step1_out_~{single_pheno}_pred.list
    >>>

    output {
        Array[File] step1_gz_files = glob("regenie_step1_out_*.gz")
        File pred_list = "regenie_step1_out_~{single_pheno}_pred.list"
        File step1_log = "regenie_step1_out_~{single_pheno}.log"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task ConcatPredLists {
    input {
        String project_id
        Array[File] pred_list_files

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 10,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    command <<<
        set -euo pipefail
        cat ~{sep=" " pred_list_files} > regenie_step1_out_pred.list
    >>>

    output {
        File merged_pred_list = "regenie_step1_out_pred.list"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task RegenieLogParserStep1 {
    input {
        String project_id
        File regenie_step1_log
        String project_name

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 10,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    command <<<
        set -euo pipefail
        RegenieLogParser.py ~{regenie_step1_log} --output ~{project_name}.step1.log
    >>>

    output {
        File parsed_log = "~{project_name}.step1.log"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

## ---------- Data Preparation Tasks ----------

task ConvertToBgen {
    input {
        String filename
        File vcf_file
        String chromosome
        String input_format  # vcf or bcf

        String? dosage_from
        String? import_dosage_certainty
        String? vcf_fixed_fid
        String? vcf_min_gq

        RuntimeAttributes runtime_attr = object {
            cpu: 4,
            memory: "8 GB",
            disk_size_gb: 200,
            max_retries: 3,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    String vcf_basename = basename(vcf_file, ".vcf.gz")
    String output_bgen = if input_format == "bcf" then basename(vcf_file, ".bcf") else vcf_basename

    command <<<
        set -euo pipefail

        DOSAGE_OPT=""
        if [ -n "~{default="" dosage_from}" ]; then
            DOSAGE_OPT="dosage=~{dosage_from}"
        fi

        DOSAGE_CERTAINTY_OPT=""
        if [ -n "~{default="" import_dosage_certainty}" ]; then
            DOSAGE_CERTAINTY_OPT="--import-dosage-certainty ~{import_dosage_certainty}"
        fi

        SAMPLE_ID_OPT="--double-id"
        if [ -n "~{default="" vcf_fixed_fid}" ]; then
            SAMPLE_ID_OPT="--const-fid ~{vcf_fixed_fid}"
        fi

        MIN_GQ_OPT=""
        if [ -n "~{default="" vcf_min_gq}" ]; then
            MIN_GQ_OPT="--vcf-min-gq ~{vcf_min_gq}"
        fi

        plink2 \
            --~{input_format} ~{vcf_file} $DOSAGE_OPT \
            --export bgen-1.2 ref-first 'bits=8' \
            $SAMPLE_ID_OPT \
            $MIN_GQ_OPT \
            $DOSAGE_CERTAINTY_OPT \
            --out ~{output_bgen} \
            --threads ~{runtime_attr.cpu} \
            --memory ~{4000}

        bgenix -g ~{output_bgen}.bgen -index
    >>>

    output {
        File bgen = "~{output_bgen}.bgen"
        File bgi = "~{output_bgen}.bgen.bgi"
        File sample = "~{output_bgen}.sample"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task MakeBgenIndex {
    input {
        String filename
        File bgen_file
        String chromosome

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 100,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    command <<<
        set -euo pipefail
        ln -sf ~{bgen_file} ~{basename(bgen_file)}
        bgenix -g ~{basename(bgen_file)} -index
    >>>

    output {
        File bgi = "~{basename(bgen_file)}.bgi"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task MakeBgenSample {
    input {
        String filename
        File bgen_file
        String chromosome

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 100,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    String bgen_basename = basename(bgen_file, ".bgen")

    command <<<
        set -euo pipefail
        qctool -g ~{bgen_file} -os ~{bgen_basename}.sample
    >>>

    output {
        File sample = "~{bgen_basename}.sample"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

## ---------- Chunking / Splitting Tasks ----------

task MakeVariantsChunks {
    input {
        String filename
        File primary_file
        File secondary_file
        File tertiary_file
        String chromosome
        File snplist_file
        File make_chunks_sql

        Int step2_gwas_chunk_size = 100000
        Array[String] chromosomes
        String snplist_type  # bgen, pgen, bed, vcf, bcf

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 50,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    command <<<
        set -euo pipefail

        ln -sf ~{snplist_file} snplist

        CHROMOSOMES_LIST="~{if chromosome == "ONE_FILE" then sep(" ", chromosomes) else chromosome}"
        USE_SQL="~{if snplist_type == "bgen" || snplist_type == "vcf" || snplist_type == "bcf" then "TRUE" else "FALSE"}"
        POS_IDX="~{if snplist_type == "pgen" then "2" else "4"}"

        if [[ "$USE_SQL" == "TRUE" ]]; then
            chromosome_sql=""
            for c in $CHROMOSOMES_LIST; do
                if [[ $chromosome_sql == "" ]]; then
                    chromosome_sql="'$c'"
                else
                    chromosome_sql="$chromosome_sql,'$c'"
                fi
            done
            sed -e 's/%CHUNK_SIZE%/~{step2_gwas_chunk_size}/' -e "s/%CHROMOSOMES%/$chromosome_sql/" ~{make_chunks_sql} > task_configured.sql
            sqlite3 snplist < task_configured.sql
            mv intervals.txt ~{filename}.GWAS-chunks.txt
        else
            grep -v "#" snplist | awk '{print $1, $'$POS_IDX' >> $1".snps"}'
            for c in $CHROMOSOMES_LIST; do
                if [ -f ${c}.snps ]; then
                    awk 'NR == 1 {start=$2}; NR > 1 && !(NR%~{step2_gwas_chunk_size}) {print $1":"start"-"$2; start = $2+1}; END { if (NR%~{step2_gwas_chunk_size}) {print $1":"start"-"$2} }' ${c}.snps > ${c}.intervals
                fi
            done
            cat *.intervals | sort -V > ~{filename}.GWAS-chunks.txt
        fi
    >>>

    output {
        File chunks_file = "~{filename}.GWAS-chunks.txt"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task MakeGenesChunks {
    input {
        String filename
        File primary_file
        File secondary_file
        File tertiary_file
        String chromosome
        File set_list_file

        Int step2_rarevar_chunk_size = 200
        Array[String] chromosomes

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 50,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    command <<<
        set -euo pipefail

        CHROMOSOMES_LIST="~{if chromosome == "ONE_FILE" then sep(" ", chromosomes) else chromosome}"

        grep -v "#" ~{set_list_file} | awk '{print $1 >> $2".genes"}'
        for c in $CHROMOSOMES_LIST; do
            if [ -f ${c}.genes ]; then
                split -l ~{step2_rarevar_chunk_size} -d ${c}.genes ~{filename}.rarevar-chunks.chr${c}.
            fi
        done
    >>>

    output {
        Array[File] chunk_files = glob("~{filename}.rarevar-chunks.chr*")
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

## ---------- Regenie Step 2 Tasks ----------

task RegenieStep2Gwas {
    input {
        String project_id
        File phenotypes_file
        PhenoMeta pheno_meta
        File covariates_file
        CovarMeta covar_meta
        AccessoryFiles accessory_files
        Array[File] step1_predictions
        String filename
        File bed_bgen_pgen
        File bim_bgi_pvar
        File fam_sample_psam
        String chromosome
        String chunk  # chr:start-end or SINGLE_CHUNK
        Array[String] chromosomes_list

        String genotypes_imputed_format
        Int regenie_bsize_step2 = 400
        String regenie_gwas_min_mac = "50"
        String regenie_min_imputation_score = "0.00"
        Boolean phenotypes_delete_missings = false
        Boolean regenie_skip_predictions = false
        Boolean regenie_ref_first_step2 = true
        Boolean regenie_firth = true
        Boolean regenie_firth_approx = true
        String regenie_range = ""
        Int maxCatLevels = 10
        String? additional_geno_format

        RuntimeAttributes runtime_attr = object {
            cpu: 8,
            memory: "16 GB",
            disk_size_gb: 200,
            max_retries: 3,
            docker: "ghcr.io/rgcgithub/regenie/regenie:v4.0.gz"
        }
    }

    String format = if genotypes_imputed_format == "vcf" || genotypes_imputed_format == "bcf" then "bgen" else genotypes_imputed_format
    String fileprefix = basename(bed_bgen_pgen, ".bgen")
    String extension = if genotypes_imputed_format == "bgen" || genotypes_imputed_format == "vcf" || genotypes_imputed_format == "bcf" then ".bgen" else ""

    command <<<
        set -euo pipefail

        # Stage step1 predictions in current directory
        for f in ~{sep=" " step1_predictions}; do
            ln -sf "$f" .
        done

        # Stage genotype files
        ln -sf ~{bed_bgen_pgen} .
        ln -sf ~{bim_bgi_pvar} .
        ln -sf ~{fam_sample_psam} .

        SPLIT_REGION=""
        if [ "~{chunk}" != "SINGLE_CHUNK" ]; then
            SPLIT_REGION="--range ~{chunk}"
        fi

        BGEN_SAMPLE=""
        if [ "~{format}" = "bgen" ]; then
            BGEN_SAMPLE="--sample ~{basename(fam_sample_psam)}"
        fi

        TEST=""
        if [ "~{pheno_meta.model}" != "additive" ]; then
            TEST="--test ~{pheno_meta.model}"
        fi

        FIRTH_APPROX=""
        if [ "~{regenie_firth_approx}" = "true" ]; then
            FIRTH_APPROX="--approx"
        fi

        FIRTH=""
        if [ "~{regenie_firth}" = "true" ]; then
            FIRTH="--firth $FIRTH_APPROX"
        fi

        BINARY_TRAIT=""
        if [ "~{pheno_meta.binary}" = "true" ]; then
            BINARY_TRAIT="--bt $FIRTH"
        fi

        RANGE=""
        if [ -n "~{regenie_range}" ]; then
            RANGE="--range ~{regenie_range}"
        fi

        EXTRACT_SNPS=""
        ~{if defined(accessory_files.extract_snps_list) then 'EXTRACT_SNPS="--extract ' + select_first([accessory_files.extract_snps_list, ""]) + '"' else ''}

        COVARIANTS=""
        if [ "~{basename(covariates_file)}" != "NO_COV_FILE" ]; then
            COVARIANTS="--covarFile ~{covariates_file} --covarColList ~{covar_meta.cols}"
        fi

        CAT_COVARIATES=""
        if [ -n "~{covar_meta.cat_cols}" ] && [ "~{covar_meta.cat_cols}" != "NA" ]; then
            CAT_COVARIATES="--catCovarList ~{covar_meta.cat_cols}"
        fi

        DELETE_MISSING_DATA=""
        if [ "~{phenotypes_delete_missings}" = "true" ]; then
            DELETE_MISSING_DATA="--strict"
        fi

        PREDICTIONS="--pred regenie_step1_out_pred.list"
        if [ "~{regenie_skip_predictions}" = "true" ]; then
            PREDICTIONS="--ignore-pred"
        fi

        REF_FIRST=""
        if [ "~{regenie_ref_first_step2}" = "true" ]; then
            REF_FIRST="--ref-first"
        fi

        CHROMOSOME=""
        if [ "~{chromosome}" != "ONE_FILE" ]; then
            CHROMOSOME="--chr ~{chromosome}"
        fi

        INTERACTION_COV=""
        if [ -n "~{default="" covar_meta.gxe}" ] && [ "~{default="" covar_meta.gxe}" != "NA" ]; then
            INTERACTION_COV="--interaction ~{covar_meta.gxe}"
        fi

        INTERACTION_SNP=""
        if [ -n "~{default="" covar_meta.gxg}" ] && [ "~{default="" covar_meta.gxg}" != "NA" ]; then
            INTERACTION_SNP="--interaction-snp ~{covar_meta.gxg}"
        fi

        regenie \
            --step 2 \
            --~{format} ~{basename(bed_bgen_pgen)} \
            --chrList ~{sep="," chromosomes_list} \
            --phenoFile ~{phenotypes_file} \
            --phenoColList ~{pheno_meta.cols} \
            --bsize ~{regenie_bsize_step2} \
            $PREDICTIONS \
            --threads ~{runtime_attr.cpu} \
            --minMAC ~{regenie_gwas_min_mac} \
            --minINFO ~{regenie_min_imputation_score} \
            --gz \
            $SPLIT_REGION \
            $BINARY_TRAIT \
            $TEST \
            $BGEN_SAMPLE \
            $CHROMOSOME \
            $RANGE \
            $EXTRACT_SNPS \
            $COVARIANTS \
            $CAT_COVARIATES \
            $DELETE_MISSING_DATA \
            $REF_FIRST \
            --maxCatLevels ~{maxCatLevels} \
            $INTERACTION_COV \
            $INTERACTION_SNP \
            --out ~{project_id}_~{filename}_~{chunk}
    >>>

    output {
        Array[File] regenie_results = glob("*.regenie.gz")
        File step2_log = "~{project_id}_~{filename}_~{chunk}.log"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task RegenieStep2Rarevars {
    input {
        String project_id
        File phenotypes_file
        PhenoMeta pheno_meta
        File covariates_file
        CovarMeta covar_meta
        AccessoryFiles accessory_files
        Array[File] step1_predictions
        String filename
        File bed_bgen_pgen
        File bim_bgi_pvar
        File fam_sample_psam
        String chromosome
        File? gene_chunk_file  # gene subset file for --extract-sets (when splitting)
        Int task_index
        Array[String] chromosomes_list

        File rarevars_set_list
        File rarevars_anno_file
        File rarevars_mask_file

        String genotypes_rarevar_format
        Int regenie_bsize_step2 = 400
        String regenie_rarevar_min_mac = "1"
        String rarevars_aaf_bins = "0.01,0.05"
        String rarevars_vc_test = "skat,skato,acatv,acato"
        String? rarevars_joint_test
        String? rarevars_vc_maxAAF
        String? regenie_build_mask
        Boolean rarevars_write_mask_snplist = false
        Boolean phenotypes_delete_missings = false
        Boolean regenie_skip_predictions = false
        Boolean regenie_ref_first_step2 = true
        Boolean regenie_firth = true
        Boolean regenie_firth_approx = true
        String regenie_range = ""
        Int maxCatLevels = 10
        String? additional_geno_format

        RuntimeAttributes runtime_attr = object {
            cpu: 8,
            memory: "42 GB",
            disk_size_gb: 200,
            max_retries: 3,
            docker: "ghcr.io/rgcgithub/regenie/regenie:v4.0.gz"
        }
    }

    String format = if genotypes_rarevar_format == "vcf" || genotypes_rarevar_format == "bcf" then "bgen" else genotypes_rarevar_format
    String fileprefix = basename(bed_bgen_pgen, ".bgen")
    String extension = if genotypes_rarevar_format == "bgen" || genotypes_rarevar_format == "vcf" || genotypes_rarevar_format == "bcf" then ".bgen" else ""

    command <<<
        set -euo pipefail

        # Stage step1 predictions in current directory
        for f in ~{sep=" " step1_predictions}; do
            ln -sf "$f" .
        done
        ln -sf ~{bed_bgen_pgen} .
        ln -sf ~{bim_bgi_pvar} .
        ln -sf ~{fam_sample_psam} .

        SPLIT_GENES=""
        if [ "~{defined(gene_chunk_file)}" = "true" ]; then
            SPLIT_GENES="--extract-sets ~{gene_chunk_file}"
        fi

        CHROMOSOME=""
        if [ "~{chromosome}" != "ONE_FILE" ]; then
            CHROMOSOME="--chr ~{chromosome}"
        fi

        BGEN_SAMPLE=""
        if [ "~{format}" = "bgen" ]; then
            BGEN_SAMPLE="--sample ~{basename(fam_sample_psam)}"
        fi

        BUILD_MASK=""
        if [ -n "~{default="" regenie_build_mask}" ]; then
            BUILD_MASK="--build-mask ~{regenie_build_mask}"
        fi

        FIRTH_APPROX=""
        if [ "~{regenie_firth_approx}" = "true" ]; then
            FIRTH_APPROX="--approx"
        fi

        FIRTH=""
        if [ "~{regenie_firth}" = "true" ]; then
            FIRTH="--firth $FIRTH_APPROX"
        fi

        BINARY_TRAIT=""
        if [ "~{pheno_meta.binary}" = "true" ]; then
            BINARY_TRAIT="--bt $FIRTH"
        fi

        COVARIANTS=""
        if [ "~{basename(covariates_file)}" != "NO_COV_FILE" ]; then
            COVARIANTS="--covarFile ~{covariates_file} --covarColList ~{covar_meta.cols}"
        fi

        CAT_COVARIATES=""
        if [ -n "~{covar_meta.cat_cols}" ] && [ "~{covar_meta.cat_cols}" != "NA" ]; then
            CAT_COVARIATES="--catCovarList ~{covar_meta.cat_cols}"
        fi

        DELETE_MISSING_DATA=""
        if [ "~{phenotypes_delete_missings}" = "true" ]; then
            DELETE_MISSING_DATA="--strict"
        fi

        PREDICTIONS="--pred regenie_step1_out_pred.list"
        if [ "~{regenie_skip_predictions}" = "true" ]; then
            PREDICTIONS="--ignore-pred"
        fi

        REF_FIRST=""
        if [ "~{regenie_ref_first_step2}" = "true" ]; then
            REF_FIRST="--ref-first"
        fi

        VC_TESTS=""
        if [ -n "~{rarevars_vc_test}" ]; then
            VC_TESTS="--vc-tests ~{rarevars_vc_test}"
        fi

        JOINT_TESTS=""
        if [ -n "~{default="" rarevars_joint_test}" ]; then
            JOINT_TESTS="--joint ~{rarevars_joint_test}"
        fi

        VC_MAXAAF=""
        if [ -n "~{default="" rarevars_vc_maxAAF}" ]; then
            VC_MAXAAF="--vc-maxAAF ~{rarevars_vc_maxAAF}"
        fi

        WRITE_MASK_SNPLIST=""
        if [ "~{rarevars_write_mask_snplist}" = "true" ]; then
            WRITE_MASK_SNPLIST="--write-mask-snplist"
        fi

        RANGE=""
        if [ -n "~{regenie_range}" ]; then
            RANGE="--range ~{regenie_range}"
        fi

        regenie \
            --step 2 \
            --~{format} ~{basename(bed_bgen_pgen)} \
            --anno-file ~{rarevars_anno_file} \
            --set-list ~{rarevars_set_list} \
            --mask-def ~{rarevars_mask_file} \
            --phenoFile ~{phenotypes_file} \
            --phenoColList ~{pheno_meta.cols} \
            --bsize ~{regenie_bsize_step2} \
            $PREDICTIONS \
            --threads ~{runtime_attr.cpu} \
            --gz \
            --aaf-bins ~{rarevars_aaf_bins} \
            --minMAC ~{regenie_rarevar_min_mac} \
            $CHROMOSOME \
            $SPLIT_GENES \
            $VC_TESTS \
            $JOINT_TESTS \
            $VC_MAXAAF \
            $BINARY_TRAIT \
            $BGEN_SAMPLE \
            $RANGE \
            $COVARIANTS \
            $CAT_COVARIATES \
            $DELETE_MISSING_DATA \
            $PREDICTIONS \
            $REF_FIRST \
            --maxCatLevels ~{maxCatLevels} \
            $BUILD_MASK \
            $WRITE_MASK_SNPLIST \
            --out ~{project_id}_~{filename}_~{task_index}
    >>>

    output {
        Array[File] regenie_results = glob("*.regenie.gz")
        File step2_log = "~{project_id}_~{filename}_~{task_index}.log"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task RegenieLogParserStep2 {
    input {
        String project_id
        Array[File] regenie_step2_logs

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 10,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    command <<<
        set -euo pipefail
        RegenieLogParser.py ~{regenie_step2_logs[0]} --output ~{project_id}.step2.log
    >>>

    output {
        File parsed_log = "~{project_id}.step2.log"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

## ---------- Results Processing Tasks ----------

task ConcatStep2Results {
    input {
        String project_id
        String phenotype
        Array[File] regenie_gz_files
        Boolean rarevar_results = false

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "8 GB",
            disk_size_gb: 100,
            max_retries: 3,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    String suffix = if rarevar_results then "rarevars" else "gwas"
    Int n_head_lines = if rarevar_results then 2 else 1

    command <<<
        set -euo pipefail
        mkdir -p tmp_sort

        for f in ~{sep=" " regenie_gz_files}; do
            zcat "$f" | tail -n+~{n_head_lines + 1} >> ~{phenotype}.tmp
        done

        headerfile="~{regenie_gz_files[0]}"
        zcat "$headerfile" | head -n ~{n_head_lines} | sed 's/ /\t/g' > header.txt

        (cat header.txt && sed 's/ /\t/g' ~{phenotype}.tmp | sort -S 6G -k1,1V -k2,2n -T tmp_sort) | bgzip -c > ~{phenotype}.~{suffix}.regenie.gz

        sort_exit=${PIPESTATUS[0]}
        if [[ $sort_exit -ne 0 ]]; then
            echo "ERROR CODE: $sort_exit"
            exit $sort_exit
        fi

        tabix -f -b 2 -e 2 -s 1 -S 1 ~{phenotype}.~{suffix}.regenie.gz
    >>>

    output {
        File merged_results = "~{phenotype}.~{suffix}.regenie.gz"
        File merged_results_tbi = "~{phenotype}.~{suffix}.regenie.gz.tbi"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task FilterResults {
    input {
        String project_id
        String phenotype
        File regenie_result_gz
        Float annotation_min_log10p = 7.3
        Boolean rarevar_results = false

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "4 GB",
            disk_size_gb: 50,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    Int n_head_lines = if rarevar_results then 2 else 1
    String grep_rarevar = if rarevar_results then "tail -n+2 | " else ""

    command <<<
        set -euo pipefail
        colnum=$(zcat ~{regenie_result_gz} | ~{grep_rarevar} head -1 | tr "\t" "\n" | cat -n | grep "LOG10P" | cut -f1 | sed -e 's/ //g')
        zcat ~{regenie_result_gz} | awk -v colnum="$colnum" 'NR <= ~{n_head_lines} {print ;}; NR > ~{n_head_lines} && $colnum >= ~{annotation_min_log10p}' > ~{basename(regenie_result_gz, ".gz")}.filtered
        gzip ~{basename(regenie_result_gz, ".gz")}.filtered
    >>>

    output {
        File filtered_results = "~{basename(regenie_result_gz, ".gz")}.filtered.gz"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task AnnotateFiltered {
    input {
        String project_id
        String phenotype
        File regenie_merged
        File genes_bed
        Int annotation_interval_kb = 25

        File collapse_closegenes_py

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "4 GB",
            disk_size_gb: 50,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    String merged_basename = basename(regenie_merged, ".gz")

    command <<<
        set -euo pipefail
        mkdir -p work

        # Save original header
        zcat ~{regenie_merged} | head -1 > header.txt

        # Sort and transform to bed file
        zcat ~{regenie_merged} | tail -n+2 | sort -T work -k1,1V -k2,2n | awk '{$2 = $2-1 OFS $2} 1' OFS='\t' > ~{merged_basename}.sorted.bed

        # Generate genome file
        cut -f1 ~{genes_bed} | uniq | awk '{OFS="\t"}; {print $1,"1"}' > genome.txt

        # Annotate closest gene with bedtools
        bedtools closest -a ~{merged_basename}.sorted.bed -b ~{genes_bed} -d -g genome.txt > ~{merged_basename}.annotated.bed
        rm ~{merged_basename}.sorted.bed

        # Generate interval around SNPs for gene annotation
        INTERVAL_BP=$((~{annotation_interval_kb} * 1000))
        awk -v interval="$INTERVAL_BP" '{OFS="\t"};{$4=$3"_%SEP%_"$4; $2=$2-interval; $3=$3+interval}; {print $0}' ~{merged_basename}.annotated.bed \
            | awk '{OFS="\t"}; $2 < 0 {$2 = 0}; {print ;}' \
            | bedtools intersect -a stdin -b ~{genes_bed} -loj \
            | cut -f1,4- | sed 's/_%SEP%_/\t/' \
            | awk '{$2 = $2-1 OFS $2} 1' OFS='\t' \
            > ~{merged_basename}.closegenes.bed

        # Merge results
        python ~{collapse_closegenes_py} ~{merged_basename}.closegenes.bed
        rm ~{merged_basename}.closegenes.bed

        # Remove duplicate column
        cut -f1,3- closegenes_collapsed.tsv > ~{merged_basename}.final.fixed.bed
        rm closegenes_collapsed.tsv

        # Write extended header
        (cat header.txt | sed " 1 s/.*/&\tCLOSEST_GENE_CHROMOSOME\tCLOSEST_GENE_START\tCLOSEST_GENE_END\tCLOSEST_GENE_NAME\tCLOSEST_GENE_DISTANCE\tGENES_~{annotation_interval_kb}KB/" && cat ~{merged_basename}.final.fixed.bed) > ~{merged_basename}.final.header.bed
        rm ~{merged_basename}.final.fixed.bed

        # Sort by p-value
        (cat ~{merged_basename}.final.header.bed | head -n 1 && cat ~{merged_basename}.final.header.bed | tail -n +2 | sort -T work -k13,13gr) | gzip > ~{merged_basename}.annotated.txt.gz
    >>>

    output {
        File annotated_results = "~{merged_basename}.annotated.txt.gz"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task ProcessRarevarResults {
    input {
        String project_id
        String phenotype
        File regenie_result_gz

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "4 GB",
            disk_size_gb: 50,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    command <<<
        set -euo pipefail
        process_regenie_rarevar.py ~{regenie_result_gz}
    >>>

    output {
        File processed_results = "~{basename(regenie_result_gz, ".gz")}.correctedP.gz"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

## ---------- Clumping Tasks ----------

task ConvertToBed {
    input {
        String chromosome
        File bgen_pgen
        File bgi_pvar
        File sample_psam
        String genotypes_imputed_format

        RuntimeAttributes runtime_attr = object {
            cpu: 4,
            memory: "8 GB",
            disk_size_gb: 200,
            max_retries: 3,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    String fileprefix = basename(bgen_pgen, ".bgen")
    String format = if genotypes_imputed_format == "pgen" then "pfile" else "bgen"
    String extension = if genotypes_imputed_format == "bgen" || genotypes_imputed_format == "vcf" || genotypes_imputed_format == "bcf" then ".bgen ref-first" else ""

    command <<<
        set -euo pipefail
        ln -sf ~{bgen_pgen} .
        ln -sf ~{bgi_pvar} .
        ln -sf ~{sample_psam} .

        BGEN_SAMPLE=""
        if [ "~{genotypes_imputed_format}" = "bgen" ]; then
            BGEN_SAMPLE="--sample ~{basename(sample_psam)}"
        fi

        plink2 \
            --~{format} ~{basename(bgen_pgen)}~{extension} \
            --make-bed \
            --memory ~{4000} \
            --threads ~{runtime_attr.cpu} \
            $BGEN_SAMPLE \
            --out ~{fileprefix}
    >>>

    output {
        File bed = "~{fileprefix}.bed"
        File bim = "~{fileprefix}.bim"
        File fam = "~{fileprefix}.fam"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task MergeBedDataset {
    input {
        Array[File] bed_files
        Array[File] bim_files
        Array[File] fam_files

        RuntimeAttributes runtime_attr = object {
            cpu: 4,
            memory: "8 GB",
            disk_size_gb: 200,
            max_retries: 3,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    command <<<
        set -euo pipefail

        # Stage all files in current directory
        for f in ~{sep=" " bed_files}; do ln -sf "$f" .; done
        for f in ~{sep=" " bim_files}; do ln -sf "$f" .; done
        for f in ~{sep=" " fam_files}; do ln -sf "$f" .; done

        # Create merge list with unique basenames
        for f in ~{sep=" " bed_files}; do
            basename "$f" .bed
        done | sort -u > files_to_merge.list

        plink \
            --merge-list files_to_merge.list \
            --make-bed \
            --memory ~{4000} \
            --threads ~{runtime_attr.cpu} \
            --out ld_panel_merged
    >>>

    output {
        File merged_bed = "ld_panel_merged.bed"
        File merged_bim = "ld_panel_merged.bim"
        File merged_fam = "ld_panel_merged.fam"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task PlinkClumping {
    input {
        String project_id
        String phenotype
        File pheno_results_gz
        String chromosome
        File bed
        File bim
        File fam
        File genes_interval
        String genotypes_build

        Float clump_p1 = 5e-8
        Float clump_p2 = 1e-4
        Int clump_kb = 250
        Int annotation_interval_kb = 25

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "4 GB",
            disk_size_gb: 50,
            max_retries: 3,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    String bfile_prefix = basename(bed, ".bed")
    String output_prefix = if chromosome != "NO_SPLIT" then chromosome else bfile_prefix

    command <<<
        set -euo pipefail

        ln -sf ~{bed} ~{bfile_prefix}.bed
        ln -sf ~{bim} ~{bfile_prefix}.bim
        ln -sf ~{fam} ~{bfile_prefix}.fam

        touch ~{output_prefix}.clumped
        touch ~{output_prefix}.clumped.ranges

        zcat ~{pheno_results_gz} | awk '{OFS="\t"}; NR == 1 {print $0, "PVAL"}; NR > 1 {print $0, 10^(-$13)}' > regenie.pval

        TARGET_CHROM=""
        if [ "~{chromosome}" != "NO_SPLIT" ]; then
            TARGET_CHROM="--chr ~{chromosome}"
        fi

        plink \
            --memory ~{2000} \
            --bfile ~{bfile_prefix} \
            $TARGET_CHROM \
            --clump regenie.pval \
            --clump-p1 ~{clump_p1} \
            --clump-p2 ~{clump_p2} \
            --clump-kb ~{clump_kb} \
            --clump-r2 0.5 \
            --clump-snp-field ID \
            --clump-field PVAL \
            --clump-range ~{genes_interval} \
            --clump-range-border ~{annotation_interval_kb} \
            --out ~{output_prefix}
    >>>

    output {
        File clumped = "~{output_prefix}.clumped"
        File clumped_ranges = "~{output_prefix}.clumped.ranges"
        File clump_log = "~{output_prefix}.log"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task MergeClumpResults {
    input {
        String project_id
        String phenotype
        Array[File] chromosome_clumps
        Array[File] chromosome_ranges

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "2 GB",
            disk_size_gb: 20,
            max_retries: 1,
            docker: "htgenomeanalysisunit/gwas-nf-pipeline:0.7"
        }
    }

    command <<<
        set -euo pipefail

        # Merge clumped files
        echo -e "CHR\tF\tSNP\tBP\tP\tTOTAL\tNSIG\tS05\tS01\tS001\tS0001\tSP2" > ~{phenotype}.toploci.tsv
        for f in ~{sep=" " chromosome_clumps}; do
            tail -n+2 "$f" | tr -s " " "\t" | sed 's/^\t//g' >> toploci.tsv
        done
        sed '/^$/d' toploci.tsv | sort -k5,5g >> ~{phenotype}.toploci.tsv

        # Merge range annotations
        echo -e "CHR\tSNP\tP\tN\tPOS\tKB\tRANGES" > ~{phenotype}.toploci.annot.tsv
        for f in ~{sep=" " chromosome_ranges}; do
            tail -n+2 "$f" | tr -s " " "\t" | sed 's/^\t//g' >> toploci.annot.tsv
        done
        sed '/^$/d' toploci.annot.tsv | sort -k3,3g >> ~{phenotype}.toploci.annot.tsv
    >>>

    output {
        File toploci = "~{phenotype}.toploci.tsv"
        File annotloci = "~{phenotype}.toploci.annot.tsv"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

## ---------- Report Tasks ----------

task ReportGwas {
    input {
        String project_id
        String phenotype
        File phenotype_file_validated
        CovarMeta covar_metadata
        File regenie_merged_results
        File annotated_tophits
        File? annotated_toploci
        File report_template
        File quarto_report_css

        String project_date
        String pipeline_version
        String manhattan_annotations = "genes"
        Float annotation_min_log10p = 7.3
        Int n_top_loci_plot = 5
        Int regional_plot_window_kb = 300
        String genotypes_build

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "32 GB",
            disk_size_gb: 50,
            max_retries: 3,
            docker: "edg1983/gwas-rmd-reports:v3.0"
        }
    }

    String results_basename = basename(regenie_merged_results, ".gz")

    command <<<
        set -euo pipefail

        TOPLOCI_ARG=""
        if [ -f "~{default="" annotated_toploci}" ]; then
            TOPLOCI_ARG="-P annotated_toploci_filename:'~{annotated_toploci}'"
        fi

        quarto render ~{report_template} \
            -P project:'~{project_id}' \
            -P date:"~{project_date}" \
            -P version:'~{pipeline_version}' \
            -P sumstat_file:'~{regenie_merged_results}' \
            -P phenotype_file:'~{phenotype_file_validated}' \
            -P phenotype:'~{phenotype}' \
            -P covariates:'~{covar_metadata.cols} - ~{covar_metadata.cat_cols}' \
            -P manhattan_annotation_type:'~{manhattan_annotations}' \
            -P annotation_min_log10p:'~{annotation_min_log10p}' \
            -P annotated_tophits_filename:'~{annotated_tophits}' \
            $TOPLOCI_ARG \
            -P max_loci:'~{n_top_loci_plot}' \
            -P regional_plot_window_kb:'~{regional_plot_window_kb}' \
            -P genome_build:'~{genotypes_build}' \
            -P gwaslab_data_dir:'/gwaslab_data/' \
            --to html

        mv ~{basename(report_template, ".qmd")}.html ~{project_id}.~{results_basename}.html
    >>>

    output {
        File report_html = "~{project_id}.~{results_basename}.html"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}

task ReportRarevar {
    input {
        String project_id
        String phenotype
        File phenotype_file_validated
        CovarMeta covar_metadata
        File regenie_merged_results
        File annotated_tophits
        File report_template
        File quarto_report_css

        String project_date
        String pipeline_version
        String genotypes_build
        Float rarevar_min_log10p = 5.0
        Float rarevar_stat_test_threshold = 1.3
        String rarevar_stat_test = "BONF_bygroup"

        RuntimeAttributes runtime_attr = object {
            cpu: 1,
            memory: "32 GB",
            disk_size_gb: 50,
            max_retries: 3,
            docker: "edg1983/gwas-rmd-reports:v3.0"
        }
    }

    String results_basename = basename(regenie_merged_results, ".gz")

    command <<<
        set -euo pipefail

        quarto render ~{report_template} \
            -P project:'~{project_id}' \
            -P date:"~{project_date}" \
            -P version:'~{pipeline_version}' \
            -P sumstat_file:'~{regenie_merged_results}' \
            -P phenotype_file:'~{phenotype_file_validated}' \
            -P phenotype:'~{phenotype}' \
            -P covariates:'~{covar_metadata.cols} - ~{covar_metadata.cat_cols}' \
            -P genome_build:'~{genotypes_build}' \
            -P tophits_min_value:'~{rarevar_min_log10p}' \
            -P sig_value_threshold:'~{rarevar_stat_test_threshold}' \
            -P significance_stat_test:'~{rarevar_stat_test}' \
            -P gwaslab_data_dir:'/gwaslab_data/' \
            --to html

        mv ~{basename(report_template, ".qmd")}.html ~{project_id}.~{results_basename}.html
    >>>

    output {
        File report_html = "~{project_id}.~{results_basename}.html"
    }

    runtime {
        cpu: runtime_attr.cpu
        memory: runtime_attr.memory
        disks: "local-disk ~{runtime_attr.disk_size_gb} HDD"
        maxRetries: runtime_attr.max_retries
        docker: runtime_attr.docker
    }
}
