process PRUNE_GENOTYPED {

  publishDir "${params.logdir}/${project_id}/logs", mode: 'copy', pattern: '*.pruned.log'
  label 'process_plink2'

  input:
    tuple val(project_id), path(genotyped_qc_bim_file), path(genotyped_qc_bed_file), path(genotyped_qc_fam_file)
  
  output:
    tuple val(project_id), path("${genotyped_qc_bim_file.baseName}.pruned.bim"), path("${genotyped_qc_bim_file.baseName}.pruned.bed"),path("${genotyped_qc_bim_file.baseName}.pruned.fam"), emit: genotypes_pruned_ch
    path "${genotyped_qc_bim_file.baseName}.pruned.log"

  script:  
  """
  # Prune, filter and convert to plink
  plink2 \
    --bfile ${genotyped_qc_bim_file.baseName} \
    --double-id --maf ${params.prune_maf} \
    --indep-pairwise ${params.prune_window_kbsize} ${params.prune_step_size} ${params.prune_r2_threshold} \
    --out ${genotyped_qc_bim_file.baseName} \
    --threads ${task.cpus} \
    --memory ${task.memory.toMega()}
  plink2 \
    --bfile ${genotyped_qc_bim_file.baseName} \
    --extract ${genotyped_qc_bim_file.baseName}.prune.in \
    --double-id \
    --make-bed \
    --out ${genotyped_qc_bim_file.baseName}.pruned \
    --threads ${task.cpus} \
    --memory ${task.memory.toMega()}
  """

}
