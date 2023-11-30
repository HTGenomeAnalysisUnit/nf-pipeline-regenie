process REPORT_GWAS {
  publishDir {"${params.outdir}/${project_id}/reports/gwas"}, mode: 'copy', pattern: '*.html'

  label 'html_report'

  input:
    tuple val(project_id), path(phenotype_file_validated), val(covar_metadata), val(phenotype), path(regenie_merged_results), path(annotated_tophits), path(annotated_toploci)
    path report_template
    path quarto_report_css
 
  output:
    path "${project_id}.${regenie_merged_results.baseName}.html"

  script:
  """
  quarto render ${report_template} \
    -P project:'${project_id}' \
    -P date:"${params.project_date}" \
    -P version:'$workflow.manifest.version' \
    -P sumstat_file:'${regenie_merged_results}' \
    -P phenotype_file:'${phenotype_file_validated}' \
    -P phenotype:'${phenotype}' \
    -P covariates:'${covar_metadata.cols} - ${covar_metadata.cat_cols}' \
    -P manhattan_annotation_type:'${params.manhattan_annotations}' \
    -P annotation_min_log10p:'${params.annotation_min_log10p}' \
    -P annotated_tophits_filename:'${annotated_tophits}' \
    -P annotated_toploci_filename:'${annotated_toploci}' \
    -P max_loci:'${params.n_top_loci_plot}' \
    -P regional_plot_window_kb:'${params.regional_plot_window_kb}' \
    -P genome_build:'${params.genotypes_build}' \
    -P gwaslab_data_dir:'/gwaslab_data/' \
    --to html
    
    mv ${report_template.baseName}.html ${project_id}.${regenie_merged_results.baseName}.html
  """
}

process REPORT_RAREVAR {
  publishDir {"${params.outdir}/${project_id}/reports/rarevar"}, mode: 'copy', pattern: '*.html'

  label 'html_report'

  input:
    tuple val(project_id), path(phenotype_file_validated), val(covar_metadata), val(phenotype), path(regenie_merged_results), path(annotated_tophits)
    path report_template
    path quarto_report_css

  output:
    path "${project_id}.${regenie_merged_results.baseName}.html"

  script:
  """
  quarto render ${report_template} \
    -P project:'${project_id}' \
    -P date:"${params.project_date}" \
    -P version:'$workflow.manifest.version' \
    -P sumstat_file:'${regenie_merged_results}' \
    -P phenotype_file:'${phenotype_file_validated}' \
    -P phenotype:'${phenotype}' \
    -P covariates:'${covar_metadata.cols} - ${covar_metadata.cat_cols}' \
    -P genome_build:'${params.genotypes_build}' \
    -P tophits_min_value:'${params.rarevar_min_log10p}' \
    -P sig_value_threshold:'${params.rarevar_stat_test_threshold}' \
    -P significance_stat_test:'${params.rarevar_stat_test}' \
    -P gwaslab_data_dir:'/gwaslab_data/' \
    --to html
    
  mv ${report_template.baseName}.html ${project_id}.${regenie_merged_results.baseName}.html
  """
}


