process SUBMIT_GWAS_RUN {
  label 'submit_gwas'
  tag "$chunk_dir"
  
  stageInMode 'copy'
  errorStrategy 'ignore'
  
  maxForks params.max_gwas_submit

  input:
    file(chunk_dir)
    file(shared_config)

  output:
    file('gwas_db.list')
    
  script:
  """
  nextflow run $projectDir/main.nf \
    -profile ${params.gwas_submit_profile} \
    -c $shared_config \
    -c $chunk_dir/gwas.conf \
    -work-dir ${workDir}/nf-wf-${workflow.sessionId}/nf-subrun-${task.index} \
    --db_folder ${params.outdir}/gwas_db/${chunk_dir}

  ls ${params.outdir}/gwas_db/*/*.bcf > gwas_db.list
  """
}