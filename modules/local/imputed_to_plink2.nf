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
def dosage_certainty = params.import_dosage_certainty ? "--import-dosage-certainty ${params.import_dosage_certainty}" : ''
"""
plink2 \
  --vcf $vcf_file dosage=${params.dosage_from} \
  --make-pgen \
  --double-id \
  ${dosage_certainty} \
  --out ${vcf_file.baseName} \
  --threads ${task.cpus} \
  --memory ${task.memory.toMega()}
"""
}

process CONVERT_TO_BGEN {

  label 'process_plink2'
  publishDir "${params.outdir}/bgen", mode: 'copy'

  input:
    path vcf_file

  output:
    tuple val("${output_prefix}"), path("${output_prefix}.bgen"), path("${output_prefix}.bgen.bgi"), path("${output_prefix}.sample"), emit: genotypes_bgen

script:
output_prefix = vcf_file.name.replaceFirst('.vcf(.gz){0,1}$','')
"""
plink2 \
  --vcf $vcf_file dosage=DS \
  --export bgen-1.2 ref-first 'bits=8' \
  --const-fid 0 \
  --out ${output_prefix} \
  --threads ${task.cpus} \
  --memory ${task.memory.toMega()}

bgenix -g ${output_prefix}.bgen -index
"""
}
