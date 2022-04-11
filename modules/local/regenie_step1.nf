workflow REGENIE_STEP1_SPLIT {
  take:
    genotype_data //tuple val(genotyped_plink_filename), path(genotyped_plink_bim_file), path(genotyped_plink_bed_file), path(genotyped_plink_fam_file)
    snplist //path snplist
    sample_id //path id
    phenos //path phenotypes_file
    covars //path covariates_file
  
  main:
    SPLITL0(genotype_data, snplist, sample_id, phenos, covars)
    jobs_ch = Channel.of(1..params.step1_n_chunks)
    step1_chunks_ch = jobs_ch.combine(SPLITL0.out)
    RUNL0(step1_chunks_ch)
    RUNL1(SPLITL0.out, RUNL0.out.collect())

  emit:
    regenie_step1_out = RUNL1.out.regenie_step1_out
    regenie_step1_out_log = RUNL1.out.regenie_step1_out_log
}

process SPLITL0 {
  input:
    tuple val(genotyped_plink_filename), file(genotyped_plink_bim_file), file(genotyped_plink_bed_file), file(genotyped_plink_fam_file)
    path snplist
    path id
    path phenotypes_file
    path covariates_file

  output:
    tuple file("regenie_step1.master"), file("regenie_step1*.snplist"), val(genotyped_plink_filename), file(genotyped_plink_bim_file), file(genotyped_plink_bed_file), file(genotyped_plink_fam_file), file(snplist), file(id), file(phenotypes_file), file(covariates_file)

  script:
  def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${params.covariates_columns}" : ''
  def deleteMissings = params.phenotypes_delete_missings  ? "--strict" : ''
  def refFirst = params.regenie_ref_first  ? "--ref-first" : ''
  """
  # qcfiles path required for keep and extract (but not actually set below)
  regenie \
    --step 1 \
    --bed ${genotyped_plink_filename} \
    --extract ${snplist} \
    --keep ${id} \
    --phenoFile ${phenotypes_file} \
    --phenoColList  ${params.phenotypes_columns} \
    $covariants \
    $deleteMissings \
    $refFirst \
    --bsize ${params.regenie_bsize_step1} \
    --split-l0 regenie_step1,${params.step1_n_chunks} \
    --out regenie_step1_splitl0
  """
}

process RUNL0 {
  label 'step1_runl0'
  
  input:
    tuple val(job_n), file(master_file), file(snp_file), val(genotyped_plink_filename), file(genotyped_plink_bim_file), file(genotyped_plink_bed_file), file(genotyped_plink_fam_file), file(snplist), file(id), file(phenotypes_file), file(covariates_file)
    
  output:
    path "${master_prefix}_job${job_n}_l0_*"

  script:
  master_prefix = master_file.getSimpleName()
  def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${params.covariates_columns}" : ''
  def deleteMissings = params.phenotypes_delete_missings  ? "--strict" : ''
  def forceStep1 = params.regenie_force_step1  ? "--force-step1" : ''
  def refFirst = params.regenie_ref_first  ? "--ref-first" : ''
  def useLoocv = params.use_loocv ? "--loocv" : ''
  """
  # qcfiles path required for keep and extract (but not actually set below)
  regenie \
    --step 1 \
    --bed ${genotyped_plink_filename} \
    --extract ${snplist} \
    --keep ${id} \
    --phenoFile ${phenotypes_file} \
    --phenoColList  ${params.phenotypes_columns} \
    $covariants \
    $deleteMissings \
    $forceStep1 \
    $refFirst \
    $useLoocv \
    --threads ${task.cpus} \
    --bsize ${params.regenie_bsize_step1} \
    ${params.phenotypes_binary_trait == true ? '--bt' : ''} \
    --run-l0 ${master_file},${job_n} \
    --out regenie_step1_run-l0_${job_n}
  """
}

process RUNL1 {
  label 'step1_runl1'
  stageInMode 'copy'
  
  publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'regenie_step1_out.log'
  if (params.save_step1_predictions) {
    publishDir "${params.outdir}/regenie_step1_preds", mode: 'copy', pattern: 'regenie_step1_out*.gz'
  }

  input:
    tuple file(master_file), file(snp_file), val(genotyped_plink_filename), file(genotyped_plink_bim_file), file(genotyped_plink_bed_file), file(genotyped_plink_fam_file), file(snplist), file(id), file(phenotypes_file), file(covariates_file)
    file(runl0_files)

  output:
    path "regenie_step1_out*", emit: regenie_step1_out
    path "regenie_step1_out.log", emit: regenie_step1_out_log

  script:
  master_prefix = master_file.getSimpleName()
  def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${params.covariates_columns}" : ''
  def deleteMissings = params.phenotypes_delete_missings  ? "--strict" : ''
  def forceStep1 = params.regenie_force_step1  ? "--force-step1" : ''
  def refFirst = params.regenie_ref_first  ? "--ref-first" : ''
  def useLoocv = params.use_loocv ? "--loocv" : ''
  """
  # qcfiles path required for keep and extract (but not actually set below)
  regenie \
    --step 1 \
    --bed ${genotyped_plink_filename} \
    --extract ${snplist} \
    --keep ${id} \
    --phenoFile ${phenotypes_file} \
    --phenoColList  ${params.phenotypes_columns} \
    $covariants \
    $deleteMissings \
    $forceStep1 \
    $refFirst \
    $useLoocv \
    --threads ${task.cpus} \
    --bsize ${params.regenie_bsize_step1} \
    ${params.phenotypes_binary_trait == true ? '--bt' : ''} \
    --run-l1 ${master_file} \
    --keep-l0 --gz --verbose \
    --out regenie_step1_out
  """
}

process REGENIE_STEP1 {
  label 'step1_single'

  publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'regenie_step1_out.log'
  if (params.save_step1_predictions) {
    publishDir "${params.outdir}/regenie_step1_preds", mode: 'copy', pattern: 'regenie_step1_out*.gz'
  }

  input:
    tuple val(genotyped_plink_filename), path(genotyped_plink_bim_file), path(genotyped_plink_bed_file), path(genotyped_plink_fam_file)
    path snplist
    path id
    path phenotypes_file
    path covariates_file

  output:
    path "regenie_step1_out*", emit: regenie_step1_out
    path "regenie_step1_out.log", emit: regenie_step1_out_log

  script:
  def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${params.covariates_columns}" : ''
  def deleteMissings = params.phenotypes_delete_missings  ? "--strict" : ''
  def forceStep1 = params.regenie_force_step1  ? "--force-step1" : ''
  def refFirst = params.regenie_ref_first  ? "--ref-first" : ''
  """
  # qcfiles path required for keep and extract (but not actually set below)
  regenie \
    --step 1 \
    --bed ${genotyped_plink_filename} \
    --extract ${snplist} \
    --keep ${id} \
    --phenoFile ${phenotypes_file} \
    --phenoColList  ${params.phenotypes_columns} \
    $covariants \
    $deleteMissings \
    $forceStep1 \
    $refFirst \
    --bsize ${params.regenie_bsize_step1} \
    ${params.phenotypes_binary_trait == true ? '--bt' : ''} \
    --lowmem \
    --gz \
    --lowmem-prefix tmp_rg \
    --threads ${task.cpus} \
    --out regenie_step1_out
  """

}
