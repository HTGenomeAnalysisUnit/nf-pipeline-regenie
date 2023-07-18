process REPORT_GWAS {
  publishDir {"${params.outdir}/${project_id}/reports/gwas"}, mode: 'copy', pattern: '*.html'

  label 'html_report'

  input:
    tuple val(project_id), path(phenotype_file_validated), path(phenotype_log), path(covariate_log), path(step1_log), path(step2_log), val(phenotype), path(regenie_merged_results), path(annotated_tophits), path(annotated_toploci)
    path report_template
    path quarto_report_css
 
  output:
    path "${project_id}.${regenie_merged_results.baseName}.html"

  script:
  """
  quarto render ${report_template} \
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
    
    mv ${report_template.baseName}.html ${project_id}.${regenie_merged_results.baseName}.html
  """
}

process REPORT_RAREVAR {
  publishDir {"${params.outdir}/${project_id}/reports/rarevar"}, mode: 'copy', pattern: '*.html'

  label 'html_report'

  input:
    tuple val(project_id), path(phenotype_file_validated), path(phenotype_log), path(covariate_log), path(step1_log), path(step2_log), val(phenotype), path(regenie_merged_results), path(annotated_tophits)
    path report_template
    path quarto_report_css

  output:
    path "${project_id}.${regenie_merged_results.baseName}.html"

  script:
      def annotation_as_string = params.manhattan_annotation_enabled.toString().toUpperCase()

  """
  quarto render ${report_template} \
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
    -P genome_build:'${params.genotypes_build}' \
    -P tophits_min_value:'${params.rarevar_min_log10p}' \
    -P sig_value_threshold:'${params.rarevar_stat_test_threshold}' \
    -P significance_stat_test:'${params.rarevar_stat_test}' \
    --to html
    
  mv ${report_template.baseName}.html ${project_id}.${regenie_merged_results.baseName}.html
  """
}


