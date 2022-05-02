process PREPARE_PROJECTS {
  label 'prepare_project'
  publishDir "${params.outdir}", mode: 'copy', pattern: 'master_table.tsv'

  input:
    file(pheno_chunker_r)
    file(prepare_projects_py)
    file(traits_table)
    file(models_table)
    file(fam_file)
    file(template_config)

  output:
    path "chunk_*", emit: chunks
    path "master_table.tsv", emit: master_table

  script:
  def chunk_size = params.pheno_chunk_size != null ? "-s ${params.pheno_chunk_size}" : ""
  def missing_tolerance = params.missing_tolerance != null ? "-t ${params.missing_tolerance}" : ""
  """
  Rscript --vanilla $pheno_chunker_r \
    --file $traits_table \
    --model_file $models_table \
    --sample_file $fam_file \
    $chunk_size \
    $missing_tolerance

  python $prepare_projects_py \
    master_table.tsv \
    $template_config \
    ${params.outdir}
  """
}