process FILTER_RESULTS {

  tag "${phenotype}"

  input:
    tuple val(phenotype), path(regenie_result_gz)

  output:
    tuple val(phenotype), path("${regenie_result_gz.baseName}.filtered.gz"), emit: results_filtered
    //tuple val(phenotype), path("${regenie_chromosomes}"), emit: results

  """
  zcat $regenie_result_gz | awk 'NR == 1 {print ;}; NR > 1 && \$13 >= ${params.annotation_min_log10p}' > ${regenie_result_gz.baseName}.filtered
  #todo: CSVWriter for gzip
  gzip ${regenie_result_gz.baseName}.filtered
  """

}
