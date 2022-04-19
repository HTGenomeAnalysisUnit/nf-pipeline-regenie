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
    tuple val("${imputed_vcf_file.baseName}"), path("${imputed_vcf_file.baseName}.bgen"), path("${imputed_vcf_file.baseName}.bgen.bgi"), emit: imputed_bgen

script:
"""
plink2 \
  --vcf $imputed_vcf_file dosage=DS \
  --export bgen-1.2 'bits=8' \
  --double-id \
  --out ${imputed_vcf_file.baseName} \
  --threads ${task.cpus} \
  --memory ${task.memory.toMega()}

bgenix -g ${imputed_vcf_file.baseName}.bgen -index
"""
}
