process SUBMIT_GWAS_RUN {
  label 'submit_gwas'
  tag "$chunk_dir"
  
  stageInMode 'copy'
  errorStrategy 'ignore'
  
  maxForks params.max_gwas_submit

  input:
    file(chunk_dir)
    file(shared_config)

  script:
  """
  nextflow run $projectDir/main.nf \
    -profile ${params.gwas_submit_profile} \
    -c $shared_config \
    -c $chunk_dir/gwas.conf \
    -work-dir ${workDir}/nf-wf-${workflow.sessionId}/nf-subrun-${task.index}
    --db_folder ${params.outdir}/gwas_db
  """
}