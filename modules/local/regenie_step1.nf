workflow REGENIE_STEP1_SPLIT {
  take:
    project_data //[project_id, pheno_file, pheno_meta(cols, binary, model), covar_file, covar_meta(cols, cat_cols), file_bim, file_bed, file_fam]
  
  main:
    SPLITL0(project_data)
    jobs_ch = Channel.of(1..params.step1_n_chunks)
    step1_chunks_ch = jobs_ch.combine(SPLITL0.out)
    RUNL0(step1_chunks_ch)

    l1_input_ch = RUNL0.out.groupTuple(size: params.step1_n_chunks)
      .join(SPLITL0.out)
    RUNL1(l1_input_ch)

  emit:
    regenie_step1_out = RUNL1.out.regenie_step1_out
    regenie_step1_out_log = RUNL1.out.regenie_step1_out_log
}

process SPLITL0 {
  input:
    tuple val(project_id), path(phenotypes_file), val(pheno_meta), path(covariates_file), val(covar_meta), path(file_bim), path(file_bed), path(file_fam)

  output:
    tuple val(project_id), path(phenotypes_file), val(pheno_meta), path(covariates_file), val(covar_meta), path("regenie_step1.master"), path("regenie_step1*.snplist"), val(genotyped_plink_filename), path(file_bim), path(file_bed), path(file_fam)

  script:
  def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${covar_meta.cols}" : ''
  def cat_covariates = covar_meta.cat_cols == '' || covar_meta.cat_cols == 'NA' ? '' : "--catCovarList ${covar_meta.cat_cols}"
  def deleteMissings = params.phenotypes_delete_missings  ? "--strict" : ''
  def refFirst = params.regenie_ref_first  ? "--ref-first" : ''
  def maxCatLevels = params.maxCatLevels ? "--maxCatLevels ${params.maxCatLevels}" : ''

  """
  # qcfiles path required for keep and extract (but not actually set below)
  regenie \
    --step 1 \
    --bed ${file_bed.baseName} \
    --phenoFile ${phenotypes_file} \
    --phenoColList  ${pheno_meta.cols} \
    $covariants \
    $cat_covariates \
    $deleteMissings \
    $refFirst \
    $maxCatLevels \
    --bsize ${params.regenie_bsize_step1} \
    --split-l0 regenie_step1,${params.step1_n_chunks} \
    --out regenie_step1_splitl0
  """
}

process RUNL0 {
  label 'step1_runl0'
  
  input:
    tuple val(job_n), val(project_id), path(phenotypes_file), val(pheno_meta), path(covariates_file), val(covar_meta), path(master_file), path(snpfile), val(genotyped_plink_filename), path(file_bim), path(file_bed), path(file_fam)
    
  output:
    tuple val(project_id), path("${master_prefix}_job${job_n}_l0_*")

  script:
  master_prefix = master_file.simpleName
  def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${covar_meta.cols}" : ''
  def cat_covariates = covar_meta.cat_cols == '' || covar_meta.cat_cols == 'NA' ? '' : "--catCovarList ${covar_meta.cat_cols}"
  def deleteMissings = params.phenotypes_delete_missings  ? "--strict" : ''
  def forceStep1 = params.regenie_force_step1  ? "--force-step1" : ''
  def refFirst = params.regenie_ref_first  ? "--ref-first" : ''
  def useLoocv = params.use_loocv ? "--loocv" : ''
  def maxCatLevels = params.maxCatLevels ? "--maxCatLevels ${params.maxCatLevels}" : ''
  def binary = pheno_meta.binary == true ? '--bt' : ''

  """
  # qcfiles path required for keep and extract (but not actually set below)
  regenie \
    --step 1 \
    --bed ${file_bed.baseName} \
    --phenoFile ${phenotypes_file} \
    --phenoColList  ${pheno_meta.cols} \
    $covariants \
    $cat_covariates \
    $deleteMissings \
    $forceStep1 \
    $refFirst \
    $useLoocv \
    $maxCatLevels \
    $binary \
    --threads ${task.cpus} \
    --bsize ${params.regenie_bsize_step1} \
    --run-l0 ${master_file},${job_n} \
    --out regenie_step1_run-l0_${job_n}
  """
}

process RUNL1 {
  label 'step1_runl1'
  stageInMode 'copy'
  
  publishDir "${params.logdir}/${project_id}/logs", mode: 'copy', pattern: 'regenie_step1_out.log'
  if (params.save_step1_predictions) {
    publishDir "${params.outdir}/${project_id}/regenie_step1_preds", mode: 'copy', pattern: 'regenie_step1_out_*'
  }

  input:
    tuple val(project_id), path(runl0_files), path(phenotypes_file), val(pheno_meta), path(covariates_file), val(covar_meta), path(master_file), path(snpfile), val(genotyped_plink_filename), path(file_bim), path(file_bed), path(file_fam)

  output:
    tuple val(project_id), path("regenie_step1_out_*"), emit: regenie_step1_out
    tuple val(project_id), path("regenie_step1_out.log"), emit: regenie_step1_out_log

  script:
  master_prefix = master_file.simpleName()
  def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${covar_meta.cols}" : ''
  def cat_covariates = covar_meta.cat_cols == '' || covar_meta.cat_cols == 'NA' ? '' : "--catCovarList ${covar_meta.cat_cols}"
  def deleteMissings = params.phenotypes_delete_missings  ? "--strict" : ''
  def forceStep1 = params.regenie_force_step1  ? "--force-step1" : ''
  def refFirst = params.regenie_ref_first  ? "--ref-first" : ''
  def useLoocv = params.use_loocv ? "--loocv" : ''
  def maxCatLevels = params.maxCatLevels ? "--maxCatLevels ${params.maxCatLevels}" : ''
  def binary = pheno_meta.binary == true ? '--bt' : ''

  """
  # qcfiles path required for keep and extract (but not actually set below)
  regenie \
    --step 1 \
    --bed ${genotyped_plink_filename} \
    --phenoFile ${phenotypes_file} \
    --phenoColList  ${params.phenotypes_columns} \
    $covariants \
    $cat_covariates \
    $deleteMissings \
    $forceStep1 \
    $refFirst \
    $useLoocv \
    $maxCatLevels \
    $binary \
    --threads ${task.cpus} \
    --bsize ${params.regenie_bsize_step1} \
    --niter ${params.niter} \
    --run-l1 ${master_file} \
    --keep-l0 --gz --verbose \
    --out regenie_step1_out

    sed -i "s|\$PWD/||g" regenie_step1_out_pred.list
  """
}