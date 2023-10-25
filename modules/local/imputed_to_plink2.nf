process CONVERT_TO_PGEN {
  label 'process_plink2'
  if (params.publish) {
    publishDir "${params.outdir}", mode: 'copy'
  }
  
  input:
    tuple val(filename), path(vcf_file), val(chrom)

  output:
    tuple val(filename), path("${vcf_file.baseName}.pgen"), path("${vcf_file.baseName}.pvar"),path("${vcf_file.baseName}.psam"), val(chrom), emit: genotypes_data
    path("${vcf_file.baseName}.pvar"), emit: variants_pvar

script:
def dosage_opt = params.dosage_from ? "dosage=${params.dosage_from}" : ''
def dosage_certainty_opt = params.import_dosage_certainty ? "--import-dosage-certainty ${params.import_dosage_certainty}" : ''
def sample_id_opt = params.vcf_fixed_fid ? "--const-fid ${params.vcf_fixed_fid}" : '--double-id'
def min_gq_opt = params.vcf_min_gq ? "--vcf-min-gq ${params.vcf_min_gq}" : ''
"""
plink2 \
  --${params.input_format} $vcf_file ${dosage_opt} \
  --make-pgen 'psam-cols=fid,sex,phenos' \
  ${sample_id_opt} \
  ${min_gq_opt} \
  ${dosage_certainty_opt} \
  --out ${vcf_file.baseName} \
  --threads ${task.cpus} \
  --memory ${task.memory.toMega()}
"""
}

process CONVERT_TO_BGEN {
  label 'process_plink2'
  if (params.publish) {
    publishDir "${params.outdir}", mode: 'copy'
  }

  input:
    tuple val(filename), path(vcf_file), val(chrom)

  output:
    tuple val(filename), path("${vcf_file.baseName}.bgen"), path("${vcf_file.baseName}.bgen.bgi"), path("${vcf_file.baseName}.sample"), val(chrom), emit: genotypes_data

script:
def dosage_opt = params.dosage_from ? "dosage=${params.dosage_from}" : ''
def dosage_certainty_opt = params.import_dosage_certainty ? "--import-dosage-certainty ${params.import_dosage_certainty}" : ''
def sample_id_opt = params.vcf_fixed_fid ? "--const-fid ${params.vcf_fixed_fid}" : '--double-id'
def min_gq_opt = params.vcf_min_gq ? "--vcf-min-gq ${params.vcf_min_gq}" : ''
"""
plink2 \
  --${params.input_format} $vcf_file ${dosage_opt} \
  --export bgen-1.2 ref-first 'bits=8' \
  ${sample_id_opt} \
  ${min_gq_opt} \
  ${dosage_certainty_opt} \
  --out ${vcf_file.baseName} \
  --threads ${task.cpus} \
  --memory ${task.memory.toMega()}

bgenix -g ${vcf_file.baseName}.bgen -index
"""
}
