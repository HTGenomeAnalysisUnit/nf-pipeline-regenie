process QC_FILTER_GENOTYPED {

  publishDir "${params.logdir}/${project_id}/logs", mode: 'copy', pattern: '*.qc.log'
  label 'process_plink2'

  input:
    tuple val(project_id), path(bed_file), path(bim_file), path(fam_file), path(phenos_tsv)

  output:
    path "${bed_file.baseName}.qc.log"
    path "${bed_file.baseName}.qc.snplist", emit: genotyped_filtered_snplist_ch
    path "${bed_file.baseName}.qc.id", emit: genotyped_filtered_id_ch
    tuple val(project_id), path("${bed_file.baseName}.qc.bim"), path("${bed_file.baseName}.qc.bed"),path("${bed_file.baseName}.qc.fam"), emit: genotyped_filtered_files_ch


  """
  grep -v "FID" $phenos_tsv | cut -f1,2 > samples.list
  
  plink2 \
    --bfile ${bed_file.baseName} \
    --keep samples.list \
    --maf ${params.qc_maf} \
    --mac ${params.qc_mac} \
    --geno ${params.qc_geno} \
    --hwe ${params.qc_hwe} \
    --mind ${params.qc_mind} \
    --write-snplist --write-samples --no-id-header \
    --out ${bed_file.baseName}.qc \
    --make-bed \
    --threads ${task.cpus} \
    --memory ${task.memory.toMega()}
  """

}
