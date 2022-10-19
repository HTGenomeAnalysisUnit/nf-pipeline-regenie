process IMPUTED_TO_PLINK2 {

  label 'process_plink2'

  input:
    path imputed_vcf_file

  output:
    tuple val("${imputed_vcf_file.baseName}"), path("${imputed_vcf_file.baseName}.pgen"), path("${imputed_vcf_file.baseName}.psam"),path("${imputed_vcf_file.baseName}.pvar"), emit: imputed_plink2

script:
"""
plink2 \
  --vcf $imputed_vcf_file dosage=DS \
  --make-pgen \
  --double-id \
  --out ${imputed_vcf_file.baseName} \
  --threads ${task.cpus} \
  --memory ${task.memory.toMega()}
"""
}

process IMPUTED_TO_BGEN {

  label 'process_plink2'
  publishDir "${params.outdir}/bgen", mode: 'copy'

  input:
    path imputed_vcf_file

  output:
    tuple val("${output_prefix}"), path("${output_prefix}.bgen"), path("${output_prefix}.bgen.bgi"), path("${output_prefix}.sample"), emit: imputed_bgen

script:
output_prefix = imputed_vcf_file.name.replaceFirst('.vcf(.gz){0,1}$','')
"""
plink2 \
  --vcf $imputed_vcf_file dosage=DS \
  --export bgen-1.2 ref-first 'bits=8' \
  --const-fid 0 \
  --out ${output_prefix} \
  --threads ${task.cpus} \
  --memory ${task.memory.toMega()}

bgenix -g ${output_prefix}.bgen -index
"""
}
