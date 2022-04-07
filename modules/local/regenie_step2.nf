process REGENIE_STEP2_BYCHR {
  if (params.save_step2_logs) {
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
  }

  label "regenie2_chr"
  tag "${plink2_pgen_file.simpleName}"

  input:
	  path(step1_out)
    tuple val(filename), path(plink_bgen_file), path(bgen_index), val(chrom)
    path phenotypes_file
    path sample_file
    path covariates_file

  output:
    tuple val(filename), path("*regenie.gz"), emit: regenie_step2_out
    path "${chrom}_${filename}.log", emit: regenie_step2_out_log

  script:
    //def format = params.genotypes_imputed_format == 'bgen' ? "--bgen" : '--pgen'
    //def extension = params.genotypes_imputed_format == 'bgen' ? ".bgen" : ''
    def bgen_sample = sample_file.name != 'NO_SAMPLE_FILE' ? "--sample $sample_file" : ''
    def test = "--test $params.regenie_test"
    def firthApprox = params.regenie_firth_approx ? "--approx" : ""
    def firth = params.regenie_firth ? "--firth $firthApprox" : ""
    def binaryTrait =  params.phenotypes_binary_trait ? "--bt $firth " : ""
    def range = params.regenie_range != '' ? "--range $params.regenie_range" : ''
    def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${params.covariates_columns}" : ''
    def deleteMissingData = params.phenotypes_delete_missings  ? "--strict" : ''
    def predictions = params.regenie_skip_predictions  ? '--ignore-pred' : ""
    def refFirst = params.regenie_ref_first  ? "--ref-first" : ''

  """
  regenie \
    --step 2 \
    --bgen ${plink_bgen_file} \
    --phenoFile ${phenotypes_file} \
    --phenoColList  ${params.phenotypes_columns} \
    --bsize ${params.regenie_bsize_step2} \
    --pred regenie_step1_out_pred.list \
    --threads ${task.cpus} \
    --minMAC ${params.regenie_min_mac} \
    --minINFO ${params.regenie_min_imputation_score} \
    --gz \
    --chr $chrom \
    $binaryTrait \
    $test \
    $bgen_sample \
    $range \
    $covariants \
    $deleteMissingData \
    $predictions \
    $refFirst \
    --out ${chrom}_${filename}
  """
}

process REGENIE_STEP2_BYCHUNK {

  if (params.save_step2_logs) {
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
  }

  label "step2_bychunk"
  tag "${plink_bgen_file.simpleName}"

  input:
	  path(step1_out)
    tuple val(filename), path(plink_bgen_file), path(bgen_index), val(chunk)
    path phenotypes_file
    path sample_file
    path covariates_file

  output:
    tuple val(filename), path("*regenie.gz"), emit: regenie_step2_out
    path "${chunk}_${filename}.log", emit: regenie_step2_out_log

  script:
    //def format = params.genotypes_imputed_format == 'bgen' ? "--bgen" : '--pgen'
    //def extension = params.genotypes_imputed_format == 'bgen' ? ".bgen" : ''
    def bgen_sample = sample_file.name != 'NO_SAMPLE_FILE' ? "--sample $sample_file" : ''
    def test = "--test $params.regenie_test"
    def firthApprox = params.regenie_firth_approx ? "--approx" : ""
    def firth = params.regenie_firth ? "--firth $firthApprox" : ""
    def binaryTrait =  params.phenotypes_binary_trait ? "--bt $firth " : ""
    def range = params.regenie_range != '' ? "--range $params.regenie_range" : ''
    def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${params.covariates_columns}" : ''
    def deleteMissingData = params.phenotypes_delete_missings  ? "--strict" : ''
    def predictions = params.regenie_skip_predictions  ? '--ignore-pred' : ""
    def refFirst = params.regenie_ref_first  ? "--ref-first" : ''

  """
  regenie \
    --step 2 \
    --bgen ${plink_bgen_file} \
    --phenoFile ${phenotypes_file} \
    --phenoColList  ${params.phenotypes_columns} \
    --bsize ${params.regenie_bsize_step2} \
    --pred regenie_step1_out_pred.list \
    --threads ${task.cpus} \
    --minMAC ${params.regenie_min_mac} \
    --minINFO ${params.regenie_min_imputation_score} \
    --gz --verbose \
    --range $chunk \
    $binaryTrait \
    $test \
    $bgen_sample \
    $range \
    $covariants \
    $deleteMissingData \
    $predictions \
    $refFirst \
    --out ${chunk}_${filename}
  """
}
