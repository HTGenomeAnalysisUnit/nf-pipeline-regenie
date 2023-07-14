process REPORT_GWAS {
  publishDir {"${params.outdir}/${project_id}/reports/gwas"}, mode: 'copy', pattern: '*.html'

  label 'html_report'

  input:
    tuple val(project_id), path(phenotype_file_validated), path(phenotype_log), path(covariate_log), path(step1_log), path(step2_log), val(phenotype), path(regenie_merged_results), path(annotated_tophits), path(annotated_toploci)
    path gwas_report_template
 
  output:
    path "${project_id}.${regenie_merged_results.baseName}.html"

  script:
  """
  quarto render ${gwas_report_template} \
    -P project:'${project_id}' \
    -P date:'${params.project_date}' \
    -P version:'$workflow.manifest.version' \
    -P regenie_merged:'${regenie_merged_results}' \
    -P phenotype_file:'${phenotype_file_validated}' \
    -P phenotype:'${phenotype_file_validated}' \
    -P covariates:'${params.covariates_columns}' \
    -P regenie_step1_log:'${step1_log}' \
    -P regenie_step2_log:'${step2_log}' \
    -P phenotype_log:'${phenotype_log}' \
    -P covariate_log:'${covariate_log}' \
    -P manhattan_annotation_type:'${params.manhattan_annotations}' \
    -P annotation_min_log10p:'${params.annotation_min_log10p}' \
    -P annotated_tophits_filename:'${annotated_tophits}z' \
    -P annotated_toploci_filename:'${annotated_toploci}' \
    -P max_loci:'${params.top_loci_n_plot}' \
    -P regional_plot_window_kb:'${params.regional_plot_window_kb}' \
    -P genome_build:'${params.genotypes_build}' \
    --to html
    
    mv ${gwas_report_template.baseName}.html ${project_id}.${regenie_merged_results.baseName}.html
  """
}

process REPORT_RAREVAR {
  publishDir {"${params.outdir}/${project_id}/reports/rarevar"}, mode: 'copy', pattern: '*.html'
  publishDir {"${params.outdir}/${project_id}/results/rarevar"}, mode: 'copy', pattern: '*.gz'

  label 'html_report'

  input:
    tuple val(project_id), path(phenotype_file_validated), path(phenotype_log), path(covariate_log), path(step1_log), path(step2_log), val(phenotype), path(regenie_merged_results), path(annotated_tophits)
    path gwas_report_template

  output:
    path "${project_id}.${regenie_merged_results.baseName}.html"
    tuple path("${regenie_merged_results.baseName}.correctedP.gz"), path("${regenie_merged_results.baseName}.correctedP.gz.tbi")

  script:
      def annotation_as_string = params.manhattan_annotation_enabled.toString().toUpperCase()

  """
  Rscript -e "require( 'rmarkdown' ); render('${gwas_report_template}',
    params = list(
      project = '${project_id}',
      date = '${params.project_date}',
      version = '$workflow.manifest.version',
      regenie_merged='${regenie_merged_results}',
      regenie_filename='${regenie_merged_results.baseName}',
      phenotype_file='${phenotype_file_validated}',
      phenotype='${phenotype}',
      covariates='${params.covariates_columns}',
      phenotype_log='${phenotype_log}',
      covariate_log='${covariate_log}',
      regenie_step1_log='${step1_log}',
      regenie_step2_log='${step2_log}',
      manhattan_annotation_enabled = $annotation_as_string,
      tophits_min_value = ${params.rarevar_tophits_min_value},
      significance_stat_test = '${params.rarevar_stat_test}'
    ),
    intermediates_dir='\$PWD',
    knit_root_dir='\$PWD',
    output_file='\$PWD/${project_id}.${regenie_merged_results.baseName}.html'
  )"

  bgzip ${regenie_merged_results.baseName}.correctedP
  tabix -f -b 2 -e 2 -s 1 -S 1 ${regenie_merged_results.baseName}.correctedP.gz
  """
}


