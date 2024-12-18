workflow REGENIE_STEP1_SPLIT {
  take:
    project_data //[project_id, pheno_file, pheno_meta(cols, binary, model), covar_file, covar_meta(cols, cat_cols), [accessory_files], file_bim, file_bed, file_fam]
  
  main:
    SPLITL0(project_data)
    jobs_ch = Channel.of(1..params.step1_n_chunks)
    step1_chunks_ch = jobs_ch.combine(SPLITL0.out)
    RUNL0(step1_chunks_ch)

    l1_input_ch = RUNL0.out.groupTuple(size: params.step1_n_chunks)
      .map { tuple(it[0], it[1].flatten()) } //We need to flatten the list of files in case there are multiple phenos in a singe project
      .join(SPLITL0.out)
      .map { tuple(it[0], it[3].cols.split(',').size(), it[3].cols.split(',').toList(), it[1], it[2], it[3], it[4], it[5], it[6], it[7], it[8], it[9], it[10], it[11])}
      .transpose(by: 2)
    RUNL1(l1_input_ch)

    runl1_pred_lists_ch = RUNL1.out.regenie_step1_out_pred_list
      .map{ project, n_phenos, step1_list_files -> [groupKey(project, n_phenos), n_phenos, step1_list_files] }
      .groupTuple()
      .map{ project, n_phenos, step1_list_files -> [project, n_phenos[0], step1_list_files]}
    
    CONCAT_PRED_LISTS(runl1_pred_lists_ch)

    output_log_ch = RUNL1.out.regenie_step1_out_log.groupTuple().map{ tuple(it[0], it[1][0]) }
    output_step1_out_ch = RUNL1.out.regenie_step1_out_gz
      .mix(CONCAT_PRED_LISTS.out)
      .map{ project, n_phenos, step1_files -> [groupKey(project, n_phenos+1), step1_files] }
      .groupTuple()
  
  emit:
    regenie_step1_out = output_step1_out_ch
    regenie_step1_out_log = output_log_ch
}

process SPLITL0 {
  label 'regenie'

  input:
    tuple val(project_id), path(phenotypes_file), val(pheno_meta), path(covariates_file), val(covar_meta), path(accessory_files), path(file_bim), path(file_bed), path(file_fam)

  output:
    tuple val(project_id), path(phenotypes_file), val(pheno_meta), path(covariates_file), val(covar_meta), path(accessory_files), path("regenie_step1.master"), path("regenie_step1*.snplist"), path(file_bim), path(file_bed), path(file_fam)

  script:
  def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${covar_meta.cols}" : ''
  //def make_no_cov_file = covariates_file.name == 'NO_COV_FILE' ? "unlink NO_COV_FILE; touch NO_COV_FILE" : ''
  def cat_covariates = !covar_meta.cat_cols || covar_meta.cat_cols == '' || covar_meta.cat_cols == 'NA' ? '' : "--catCovarList ${covar_meta.cat_cols}"
  def deleteMissings = params.phenotypes_delete_missings  ? "--strict" : ''
  def refFirst = params.regenie_ref_first  ? "--ref-first" : ''
  def maxCatLevels = params.maxCatLevels ? "--maxCatLevels ${params.maxCatLevels}" : ''
  def condition_list = accessory_files[0].name != 'NO_CONDITION_FILE' ? "--condition-list ${accessory_files[0]}" : ''
  def additional_geno_fileprefix = accessory_files[1].baseName
  def additional_geno_extension = params.additional_geno_format == 'bgen' ? '.bgen' : ''
  def additional_genotypes = accessory_files[1].name != 'NO_ADDITIONAL_GENO_FILE' && accessory_files[0].name != 'NO_CONDITION_FILE' ? "--condition-file ${params.additional_geno_format},${additional_geno_fileprefix}${additional_geno_extension}" : ''
  def additional_sample_file = accessory_files[1].name != 'NO_ADDITIONAL_GENO_FILE' && accessory_files[0].name != 'NO_CONDITION_FILE' && params.additional_geno_format == 'bgen' ? "--condition-file-sample ${accessory_files[2]}" : ''

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
    $condition_list \
    $additional_genotypes \
    $additional_sample_file \
    --bsize ${params.regenie_bsize_step1} \
    --split-l0 regenie_step1,${params.step1_n_chunks} \
    --out regenie_step1_splitl0
  """
}

process RUNL0 {
  label 'regenie'
  label 'step1_runl0'
  
  input:
    tuple val(job_n), val(project_id), path(phenotypes_file), val(pheno_meta), path(covariates_file), val(covar_meta), file(accessory_files), path(master_file), path(snpfile), path(file_bim), path(file_bed), path(file_fam)
    
  output:
    tuple val(project_id), path("${master_prefix}_job${job_n}_l0_*")

  script:
  master_prefix = master_file.simpleName
  def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${covar_meta.cols}" : ''
  def cat_covariates = !covar_meta.cat_cols || covar_meta.cat_cols == '' || covar_meta.cat_cols == 'NA' ? '' : "--catCovarList ${covar_meta.cat_cols}"
  def deleteMissings = params.phenotypes_delete_missings  ? "--strict" : ''
  def forceStep1 = params.regenie_force_step1  ? "--force-step1" : ''
  def refFirst = params.regenie_ref_first  ? "--ref-first" : ''
  def useLoocv = params.use_loocv ? "--loocv" : ''
  def maxCatLevels = params.maxCatLevels ? "--maxCatLevels ${params.maxCatLevels}" : ''
  def binary = pheno_meta.binary == 'true' ? '--bt' : ''
  def condition_list = accessory_files[0].name != 'NO_CONDITION_FILE' ? "--condition-list ${accessory_files[0]}" : ''
  def additional_geno_fileprefix = accessory_files[1].baseName
  def additional_geno_extension = params.additional_geno_format == 'bgen' ? '.bgen' : ''
  def additional_genotypes = accessory_files[1].name != 'NO_ADDITIONAL_GENO_FILE' && accessory_files[0].name != 'NO_CONDITION_FILE' ? "--condition-file ${params.additional_geno_format},${additional_geno_fileprefix}${additional_geno_extension}" : ''
  def additional_sample_file = accessory_files[1].name != 'NO_ADDITIONAL_GENO_FILE' && accessory_files[0].name != 'NO_CONDITION_FILE' && params.additional_geno_format == 'bgen' ? "--condition-file-sample ${accessory_files[2]}" : ''

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
    $condition_list \
    $additional_genotypes \
    $additional_sample_file \
    --threads ${task.cpus} \
    --bsize ${params.regenie_bsize_step1} \
    --run-l0 ${master_file},${job_n} \
    --out regenie_step1_run-l0_${job_n}
  """
}

process RUNL1 {
  label 'regenie'
  label 'step1_runl1'
  //stageInMode 'copy'
  
  publishDir {"${params.logdir}/${project_id}/logs"}, mode: 'copy', pattern: 'regenie_step1_out_*.log'
  if (params.save_step1_predictions) {
    publishDir {"${params.outdir}/${project_id}/regenie_step1_preds"}, mode: 'copy', pattern: 'regenie_step1_out_*.gz'
  }

  input:
    tuple val(project_id), val(n_phenos), val(single_pheno), path(runl0_files), path(phenotypes_file), val(pheno_meta), path(covariates_file), val(covar_meta), file(accessory_files), path(master_file), path(snpfile), path(file_bim), path(file_bed), path(file_fam)

  output:
    tuple val(project_id), val(n_phenos), path("regenie_step1_out_*.gz"), emit: regenie_step1_out_gz
    tuple val(project_id), val(n_phenos), path("regenie_step1_out_${single_pheno}_pred.list"), emit: regenie_step1_out_pred_list
    tuple val(project_id), path("regenie_step1_out_${single_pheno}.log"), emit: regenie_step1_out_log

  script:
  master_prefix = master_file.simpleName
  def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${covar_meta.cols}" : ''
  def cat_covariates = !covar_meta.cat_cols || covar_meta.cat_cols == '' || covar_meta.cat_cols == 'NA' ? '' : "--catCovarList ${covar_meta.cat_cols}"
  def deleteMissings = params.phenotypes_delete_missings  ? "--strict" : ''
  def forceStep1 = params.regenie_force_step1  ? "--force-step1" : ''
  def refFirst = params.regenie_ref_first  ? "--ref-first" : ''
  def useLoocv = params.use_loocv ? "--loocv" : ''
  def maxCatLevels = params.maxCatLevels ? "--maxCatLevels ${params.maxCatLevels}" : ''
  def binary = pheno_meta.binary == 'true' ? '--bt' : ''
  def condition_list = accessory_files[0].name != 'NO_CONDITION_FILE' ? "--condition-list ${accessory_files[0]}" : ''
  def additional_geno_fileprefix = accessory_files[1].baseName
  def additional_geno_extension = params.additional_geno_format == 'bgen' ? '.bgen' : ''
  def additional_genotypes = accessory_files[1].name != 'NO_ADDITIONAL_GENO_FILE' && accessory_files[0].name != 'NO_CONDITION_FILE' ? "--condition-file ${params.additional_geno_format},${additional_geno_fileprefix}${additional_geno_extension}" : ''
  def additional_sample_file = accessory_files[1].name != 'NO_ADDITIONAL_GENO_FILE' && accessory_files[0].name != 'NO_CONDITION_FILE' && params.additional_geno_format == 'bgen' ? "--condition-file-sample ${accessory_files[2]}" : ''

  """
  # qcfiles path required for keep and extract (but not actually set below)
  regenie \
    --step 1 \
    --bed ${file_bed.baseName} \
    --phenoFile ${phenotypes_file} \
    --phenoColList ${pheno_meta.cols} \
    --l1-phenoList ${single_pheno} \
    $covariants \
    $cat_covariates \
    $deleteMissings \
    $forceStep1 \
    $refFirst \
    $useLoocv \
    $maxCatLevels \
    $binary \
    $condition_list \
    $additional_genotypes \
    $additional_sample_file \
    --threads ${task.cpus} \
    --bsize ${params.regenie_bsize_step1} \
    --niter ${params.niter} \
    --run-l1 ${master_file} \
    --keep-l0 --gz --verbose \
    --out regenie_step1_out_${single_pheno}

    sed -i "s|\$PWD/||g" regenie_step1_out_${single_pheno}_pred.list
  """
}

process CONCAT_PRED_LISTS {
  label 'small_task'

  if (params.save_step1_predictions) {
    publishDir {"${params.outdir}/${project_id}/regenie_step1_preds"}, mode: 'copy', pattern: 'regenie_step1_out_pred.list'
  }

  input:
    tuple val(project_id), val(n_phenos), path(pred_list_files)

  output:
    tuple val(project_id), val(n_phenos), path("regenie_step1_out_pred.list")
  
  script:
  """
  cat regenie_step1_out_*pred.list > final.list
  mv final.list regenie_step1_out_pred.list
  """
}
