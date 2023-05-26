process FILTER_RESULTS {
  tag "${phenotype}"

  input:
    tuple val(phenotype), path(regenie_result_gz)

  output:
    tuple val(phenotype), path("${regenie_result_gz.baseName}.filtered.gz"), emit: results_filtered
    //tuple val(phenotype), path("${regenie_chromosomes}"), emit: results

  shell:
  def n_head_lines = params.rarevar_results ? 2 : 1 
  def grep_rarevar = params.rarevar_results ? 'tail -n+2 | ' : ''
  '''
  colnum=$(zcat !{regenie_result_gz} | !{grep_rarevar} head -1 | tr " " "\\n" | cat -n | grep "LOG10P" | cut -f1 | sed -e 's/ //g')
  zcat !{regenie_result_gz} | awk "NR <= !{n_head_lines} {print ;}; NR > !{n_head_lines} && \$$colnum >= !{params.annotation_min_log10p}" > !{regenie_result_gz.baseName}.filtered
  #todo: CSVWriter for gzip
  gzip !{regenie_result_gz.baseName}.filtered
  '''
}
