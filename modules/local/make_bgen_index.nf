process MAKE_BGEN_INDEX {
  label 'process_bgenix'
  if (params.publish) {
    publishDir "${params.outdir}/bgen", mode: 'copy', pattern: '*.bgi'
  }

  input:
    path imputed_bgen_file

  output:
    tuple path("${imputed_bgen_file}"), path("${imputed_bgen_file}.bgi")

script:
"""
bgenix -g ${imputed_bgen_file} -index
"""
}