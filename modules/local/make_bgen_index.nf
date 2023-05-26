process MAKE_BGEN_INDEX {
  label 'process_bgenix'
  if (params.publish) {
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.bgi'
  }

  input:
    tuple val(filename), file(bgen_file), val(chrom)

  output:
    tuple val(filename), path(bgen_file), path("${bgen_file}.bgi"), val(chrom)

script:
"""
bgenix -g ${bgen_file} -index
"""
}