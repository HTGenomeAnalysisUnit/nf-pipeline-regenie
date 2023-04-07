process REGENIE_STEP2 {
  if (params.save_step2_logs) {
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
  }

  label "regenie2_chr"
  tag "${plink_bgen_file.simpleName}"

  input:
	  path(step1_out)
    tuple val(filename), path(plink_bgen_file), path(bgen_index), path(sample_file), val(chunk)
    path phenotypes_file
    path covariates_file

  output:
    tuple val(filename), val(chunk), path("*regenie.gz"), emit: regenie_step2_out
    path "${filename}_${chunk}.log", emit: regenie_step2_out_log

  script:
    //def format = params.genotypes_imputed_format == 'bgen' ? "--bgen" : '--pgen'
    //def extension = params.genotypes_imputed_format == 'bgen' ? ".bgen" : ''
    def split_region = params.step2_split_by == 'chunk' ? "--range $chunk" : ''
    def split_region = params.step2_split_by == 'chr' ? "--chr $chunk" : ''
    def bgen_sample = sample_file.exists() ? "--sample $sample_file" : ''
    def test = "--test $params.regenie_test"
    def firthApprox = params.regenie_firth_approx ? "--approx" : ""
    def firth = params.regenie_firth ? "--firth $firthApprox" : ""
    def binaryTrait =  params.phenotypes_binary_trait ? "--bt $firth " : ""
    def range = params.regenie_range != '' ? "--range $params.regenie_range" : ''
    def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${params.covariates_columns}" : ''
    def cat_covariates = params.covariates_cat_columns == '' || params.covariates_cat_columns == 'NA' ? '' : "--catCovarList ${params.covariates_cat_columns}"
    def deleteMissingData = params.phenotypes_delete_missings ? "--strict" : ''
    def predictions = params.regenie_skip_predictions ? '--ignore-pred' : ""
    def refFirst = params.regenie_ref_first_step2 ? "--ref-first" : ''
    def maxCatLevels = params.maxCatLevels ? "--maxCatLevels ${params.maxCatLevels}" : ''

  """
  regenie \
    --step 2 \
    --bgen ${plink_bgen_file} \
    --phenoFile ${phenotypes_file} \
    --phenoColList ${params.phenotypes_columns} \
    --bsize ${params.regenie_bsize_step2} \
    --pred regenie_step1_out_pred.list \
    --threads ${task.cpus} \
    --minMAC ${params.regenie_min_mac} \
    --minINFO ${params.regenie_min_imputation_score} \
    --gz \
    $split_region \
    $binaryTrait \
    $test \
    $bgen_sample \
    $range \
    $covariants \
    $cat_covariates \
    $deleteMissingData \
    $predictions \
    $refFirst \
    $maxCatLevels \
    --out ${filename}_${chunk}
  """
}

process REGENIE_STEP2_RAREVARS {
  if (params.save_step2_logs) {
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
  }

  label "regenie2_chr"
  tag "${plink_bgen_file.simpleName}"

  input:
	  path(step1_out)
    tuple val(filename), file(bed_or_bgen_file), file(bim_or_bgi_file), file(fam_or_sample_file), val(gene_name)
    path phenotypes_file
    path covariates_file
    path rarevars_anno_file
    path rarevars_set_list
    path rarevars_mask_file

  output:
    tuple val(filename), val(chunk), path("*regenie.gz"), emit: regenie_step2_out
    path "${filename}_${chunk}.log", emit: regenie_step2_out_log

  script:
    def format = params.sequencing_data_format in ['vcf', 'bgen'] ? "--bgen" : '--bed'
    def extension = params.sequencing_data_format in ['vcf', 'bgen'] ? '.bgen' : ''
    def bgen_sample = params.sequencing_data_format in ['vcf', 'bgen'] && fam_or_sample_file.exists() ? "--sample $fam_or_sample_file" : ''
    def test = "--test $params.regenie_test"
    def firthApprox = params.regenie_firth_approx ? "--approx" : ""
    def firth = params.regenie_firth ? "--firth $firthApprox" : ""
    def binaryTrait =  params.phenotypes_binary_trait ? "--bt $firth " : ""
    def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${params.covariates_columns}" : ''
    def cat_covariates = params.covariates_cat_columns == '' || params.covariates_cat_columns == 'NA' ? '' : "--catCovarList ${params.covariates_cat_columns}"
    def deleteMissingData = params.phenotypes_delete_missings ? "--strict" : ''
    def predictions = params.regenie_skip_predictions ? '--ignore-pred' : ""
    def refFirst = params.regenie_ref_first_step2 ? "--ref-first" : ''
    def maxCatLevels = params.maxCatLevels ? "--maxCatLevels ${params.maxCatLevels}" : ''
    def vc_tests = params.rarevars_vc_test ? "--vc-tests ${params.rarevars_vc_test}" : ''
    def vc_maxAAF = params.rarevars_vc_test ? "--vc-maxAAF ${params.rarevars_vc_maxAAF}" : ''
    def write_mask_snplist = params.rarevars_write_mask_snplist ? "--write-mask-snplist" : ''
  """
  regenie \
    --step 2 \
    ${format} ${filename}${extension} \
    --anno-file $rarevars_anno_file \
    --set-list $rarevars_set_list \
    --mask-def $rarevars_mask_file \
    --phenoFile ${phenotypes_file} \
    --phenoColList ${params.phenotypes_columns} \
    --bsize ${params.regenie_bsize_step2} \
    --pred regenie_step1_out_pred.list \
    --threads ${task.cpus} \
    --gz \
    --aaf-bins ${params.rarevars_aaf_bins} \
    $vc_tests \
    $vc_maxAAF \
    $binaryTrait \
    $test \
    $bgen_sample \
    $range \
    $covariants \
    $cat_covariates \
    $deleteMissingData \
    $predictions \
    $refFirst \
    $maxCatLevels \
    $write_mask_snplist \
    --out ${filename}_${gene_name}
  """
}


/*
- read from bed / bim / fam, bgen or VCF
- additional inputs: mask file, annotation file, set list file, 
- optional inputs: external aaf file, 
- separate minMAC for rare
- options: --aaf-bins, --vs-test, --vc-maxAAF
- optional output mask snp list: --write-mask-snplist

- split by chr based on file glob pattern
- split by gene --> with a single file --extract-setlist to extract a specific gene in the set list file
                --> with multiple files by chr we must join by chr first
- use 2 separate input dataset for common step2 and rare step2
- split present workflow into sub-workflows 
        --> processing of pheno and covars remain the same
        --> a common step1 sub-workflow 
        --> common vars step2 sub-wf 
        --> rare vars step2 sub-wf

- eventually collect results by chr, or by gene 
- prepare rare vars report 

*/