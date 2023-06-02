process REGENIE_STEP2_GWAS {
  if (params.save_step2_logs) {
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.log'
  }

  label "step2_gwas" //FIXME Need to use a params to select based on running by chunk or whole dataset
  tag "${filename}_${chrom}_${chunk}"

  input:
	  path(step1_out)
    tuple val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom), val(chunk)
    path phenotypes_file
    path covariates_file

  output:
    tuple val(filename), val(chrom), val(chunk), path("*regenie.gz"), emit: regenie_step2_out
    path "${filename}_${chrom}_${chunk}.log", emit: regenie_step2_out_log

  script:
    def format = params.genotypes_imputed_format == 'vcf' ? 'pgen' : "${params.genotypes_imputed_format}"
    def fileprefix = bed_bgen_pgen.simpleName
    def extension = params.genotypes_imputed_format == 'bgen' ? '.bgen' : ''
    def split_region = chunk == 'SINGLE_CHUNK' ? '' : "--range $chunk"
    def bgen_sample = params.genotypes_imputed_format == 'bgen' ? "--sample $fam_sample_psam" : ''
    def test = params.regenie_gwas_test != 'additive' ? "--test ${params.regenie_gwas_test}" : ''
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
    def chromosome = chrom == "ONE_FILE" ? '' : "--chr $chrom"

  """
  regenie \
    --step 2 \
    --${format} ${fileprefix}${extension} \
    --chrList ${params.chromosomes.join(',')} \
    --phenoFile ${phenotypes_file} \
    --phenoColList ${params.phenotypes_columns} \
    --bsize ${params.regenie_bsize_step2} \
    --pred regenie_step1_out_pred.list \
    --threads ${task.cpus} \
    --minMAC ${params.regenie_gwas_min_mac} \
    --minINFO ${params.regenie_min_imputation_score} \
    --gz \
    $split_region \
    $binaryTrait \
    $test \
    $bgen_sample \
    $chromosome \
    $range \
    $covariants \
    $cat_covariates \
    $deleteMissingData \
    $predictions \
    $refFirst \
    $maxCatLevels \
    --out ${filename}_${chrom}_${chunk}
  """
}

process REGENIE_STEP2_RAREVARS {
  if (params.save_step2_logs) {
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
  }

  label "step2_rarevar" //FIXME Need to use a params to select based on running by chunk or whole dataset
  tag "${filename}_${chrom}_${gene}"

  input:
	  path(step1_out)
    tuple val(filename), file(bed_bgen_pgen), file(bim_bgi_pvam), file(fam_sample_psam), val(chrom), val(gene)
    path phenotypes_file
    path covariates_file
    path rarevars_anno_file
    path rarevars_set_list
    path rarevars_mask_file

  output:
    tuple val(filename), val(chrom), val(gene), path("*regenie.gz"), emit: regenie_step2_out
    path "${filename}_${gene}.log", emit: regenie_step2_out_log

  script:
    def format = params.genotypes_rarevar_format == 'vcf' ? 'pgen' : "${params.genotypes_rarevar_format}"
    def fileprefix = bed_bgen_pgen.baseName
    def extension = params.genotypes_rarevar_format == 'bgen' ? '.bgen' : ''
    def split_genes = gene == 'SINGLE_CHUNK' ? '' : "--extract-sets $gene"
    def chromosome = chrom == "ONE_FILE" ? '' : "--chr $chrom"
    def bgen_sample = params.genotypes_rarevar_format == 'bgen' ? "--sample $fam_sample_psam" : ''
    def build_mask = params.regenie_build_mask ? "--build-mask ${params.regenie_build_mask}" : ''
    def firthApprox = params.regenie_firth_approx ? "--approx" : ""
    def firth = params.regenie_firth ? "--firth $firthApprox" : ""
    def binaryTrait =  params.phenotypes_binary_trait ? "--bt $firth " : ""
    def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${params.covariates_columns}" : ''
    def cat_covariates = params.covariates_cat_columns == '' || params.covariates_cat_columns == 'NA' ? '' : "--catCovarList ${params.covariates_cat_columns}"
    def deleteMissingData = params.phenotypes_delete_missings ? "--strict" : ''
    def predictions = params.regenie_skip_predictions ? '--ignore-pred' : ""
    def refFirst = params.regenie_ref_first_step2 ? "--ref-first" : ''
    def maxCatLevels = params.maxCatLevels ? "--maxCatLevels ${params.maxCatLevels}" : ''
    def vc_tests = params.rarevars_vc_test ? "--vc-tests ${params.rarevars_vc_test.toLowerCase()}" : ''
    def vc_maxAAF = params.rarevars_vc_maxAAF ? "--vc-maxAAF ${params.rarevars_vc_maxAAF}" : ''
    def write_mask_snplist = params.rarevars_write_mask_snplist ? "--write-mask-snplist" : ''
    def range = params.regenie_range != '' ? "--range $params.regenie_range" : ''

  """
  regenie \
    --step 2 \
    --${format} ${fileprefix}${extension} \
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
    --minMAC ${params.regenie_rarevar_min_mac} \
    $chromosome \
    $split_genes \
    $vc_tests \
    $vc_maxAAF \
    $binaryTrait \
    $bgen_sample \
    $range \
    $covariants \
    $cat_covariates \
    $deleteMissingData \
    $predictions \
    $refFirst \
    $maxCatLevels \
    $build_mask \
    $write_mask_snplist \
    --out ${filename}_${chrom}_${task.index}
  """
}

