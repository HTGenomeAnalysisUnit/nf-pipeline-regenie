process FILTER_RESULTS {
  tag "${project_id}_${phenotype}"
  
  if (params.publish) {
    publishDir "${params.outdir}/${project_id}/results/${params.rarevar_results ? 'rarevar' : 'gwas'}/tophits", mode: 'copy'
  }

  input:
    tuple val(project_id), val(phenotype), path(regenie_result_gz)

  output:
    tuple val(project_id), val(phenotype), path("${regenie_result_gz.baseName}.filtered.gz"), emit: results_filtered
    //tuple val(phenotype), path("${regenie_chromosomes}"), emit: results

  shell:
  n_head_lines = params.rarevar_results ? 2 : 1 
  grep_rarevar = params.rarevar_results ? 'tail -n+2 | ' : ''
  '''
  colnum=$(zcat !{regenie_result_gz} | !{grep_rarevar} head -1 | tr "\\t" "\\n" | cat -n | grep "LOG10P" | cut -f1 | sed -e 's/ //g')
  zcat !{regenie_result_gz} | awk -v colnum="$colnum" 'NR <= !{n_head_lines} {print ;}; NR > !{n_head_lines} && $colnum >= !{params.annotation_min_log10p}' > !{regenie_result_gz.baseName}.filtered
  gzip !{regenie_result_gz.baseName}.filtered
  '''
}
