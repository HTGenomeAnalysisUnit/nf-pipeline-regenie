process REPORT_GWAS {
  publishDir "${params.outdir}/reports", mode: 'copy', pattern: '*.html'

  label 'html_report'

  input:
    tuple val(phenotype), path(regenie_merged), path(annotated_tophits), path(annotated_toploci)
    path phenotype_file_validated
    path gwas_report_template
    path phenotype_log
    path covariate_log
    path step1_log
    path step2_log

  output:
    path "${params.project}.${regenie_merged.baseName}.html"

  script:
      def annotation_as_string = params.manhattan_annotation_enabled.toString().toUpperCase()

  """
  Rscript -e "require( 'rmarkdown' ); render('${gwas_report_template}',
    params = list(
      project = '${params.project}',
      date = '${params.project_date}',
      version = '$workflow.manifest.version',
      regenie_merged='${regenie_merged}',
      regenie_filename='${regenie_merged.baseName}',
      phenotype_file='${phenotype_file_validated}',
      phenotype='${phenotype}',
      covariates='${params.covariates_columns}',
      phenotype_log='${phenotype_log}',
      covariate_log='${covariate_log}',
      regenie_step1_log='${step1_log}',
      regenie_step2_log='${step2_log}',
      plot_ylimit=${params.plot_ylimit},
      annotated_tophits_filename='${annotated_tophits}',
      annotated_toploci_filename='${annotated_toploci}',
      manhattan_annotation_enabled = $annotation_as_string,
      annotation_min_log10p = ${params.annotation_min_log10p}
    ),
    intermediates_dir='\$PWD',
    knit_root_dir='\$PWD',
    output_file='\$PWD/${params.project}.${regenie_merged.baseName}.html'
  )"
  """
}

process REPORT_RAREVAR {
  publishDir "${params.outdir}/reports", mode: 'copy', pattern: '*.html'
  publishDir "${params.outdir}/results", mode: 'copy', pattern: '*.gz'

  label 'html_report'

  input:
    tuple val(phenotype), path(regenie_merged), path(annotated_tophits)
    path phenotype_file_validated
    path gwas_report_template
    path phenotype_log
    path covariate_log
    path step1_log
    path step2_log

  output:
    path "${params.project}.${regenie_merged.baseName}.html"
    tuple path("${regenie_merged.baseName}.correctedP.gz"), path("${regenie_merged.baseName}.correctedP.gz.tbi")

  script:
      def annotation_as_string = params.manhattan_annotation_enabled.toString().toUpperCase()

  """
  Rscript -e "require( 'rmarkdown' ); render('${gwas_report_template}',
    params = list(
      project = '${params.project}',
      date = '${params.project_date}',
      version = '$workflow.manifest.version',
      regenie_merged='${regenie_merged}',
      regenie_filename='${regenie_merged.baseName}',
      phenotype_file='${phenotype_file_validated}',
      phenotype='${phenotype}',
      covariates='${params.covariates_columns}',
      phenotype_log='${phenotype_log}',
      covariate_log='${covariate_log}',
      regenie_step1_log='${step1_log}',
      regenie_step2_log='${step2_log}',
      plot_ylimit=${params.plot_ylimit},
      manhattan_annotation_enabled = $annotation_as_string,
      tophits_min_value = ${params.rarevar_tophits_min_value},
      significance_stat_test = ${params.rarevar_stat_test}
    ),
    intermediates_dir='\$PWD',
    knit_root_dir='\$PWD',
    output_file='\$PWD/${params.project}.${regenie_merged.baseName}.html'
  )"

  bgzip ${regenie_merged.baseName}.correctedP
  tabix -f -b 2 -e 2 -s 1 -S 1 ${regenie_merged.baseName}.correctedP.gz
  """
}


