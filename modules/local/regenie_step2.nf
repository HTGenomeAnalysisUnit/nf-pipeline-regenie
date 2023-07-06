process REGENIE_STEP2_GWAS {
  if (params.save_step2_logs) {
    publishDir "${params.outdir}/${project_id}/logs/step2_gwas_logs", mode: 'copy', pattern: '*.log'
  }

  label "step2_gwas" //FIXME Need to use a params to select based on running by chunk or whole dataset
  tag "${project_id}_${chrom}_${chunk}"

  input:
	  tuple val(project_id), path(phenotypes_file), val(pheno_meta), path(covariates_file), val(covar_meta), path(step1_predictions), val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom), val(chunk)

  output:
    tuple val(project_id), val(chrom), val(chunk), path("*regenie.gz"), emit: regenie_step2_out
    tuple val(project_id), path("${project_id}_${chrom}_${chunk}.log"), emit: regenie_step2_out_log

  script:
    def format = params.genotypes_imputed_format == 'vcf' ? 'pgen' : "${params.genotypes_imputed_format}"
    def fileprefix = bed_bgen_pgen.baseName
    def extension = params.genotypes_imputed_format == 'bgen' ? '.bgen' : ''
    def split_region = chunk == 'SINGLE_CHUNK' ? '' : "--range $chunk"
    def bgen_sample = params.genotypes_imputed_format == 'bgen' ? "--sample $fam_sample_psam" : ''
    def test = pheno_meta.model != 'additive' ? "--test ${pheno_meta.model}" : ''
    def firthApprox = params.regenie_firth_approx ? "--approx" : ""
    def firth = params.regenie_firth ? "--firth $firthApprox" : ""
    def binaryTrait =  pheno_meta.binary ? "--bt $firth " : ""
    def range = params.regenie_range != '' ? "--range $params.regenie_range" : ''
    def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${covar_meta.cols}" : ''
    def cat_covariates = covar_meta.cat_cols == '' || covar_meta.cat_cols == 'NA' ? '' : "--catCovarList ${covar_meta.cat_cols}"
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
    --phenoColList ${pheno_meta.cols} \
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
    --out ${project_id}_${chrom}_${chunk}
  """
}

process REGENIE_STEP2_RAREVARS {
  if (params.save_step2_logs) {
    publishDir "${params.outdir}/${project_id}/logs/step2_rarevar_logs", mode: 'copy', pattern: '*.log'
  }

  label "step2_rarevar" //FIXME Need to use a params to select based on running by chunk or whole dataset
  tag "${project_id}_chr${chrom}_chunk${task.index}"

  input:
	  tuple val(project_id), path(phenotypes_file), val(pheno_meta), path(covariates_file), covar_meta, path(step1_predictions), val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom), val(gene)

  output:
    tuple val(project_id), val(chrom), val(task.index), path("*regenie.gz"), emit: regenie_step2_out
    tuple val(project_id), path("${project_id}_${chrom}_${task.index}.log"), emit: regenie_step2_out_log

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
    def binaryTrait =  pheno_meta.binary ? "--bt $firth " : ""
    def covariants = covariates_file.name != 'NO_COV_FILE' ? "--covarFile $covariates_file --covarColList ${covar_meta.cols}" : ''
    def cat_covariates = covar_meta.cat_cols == '' || covar_meta.cat_cols == 'NA' ? '' : "--catCovarList ${covar_meta.cat_cols}"
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
    --phenoColList ${pheno_meta.cols} \
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
    --out ${project_id}_${chrom}_${task.index}
  """
}

