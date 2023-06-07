process MAKE_BGEN_SAMPLE {
  label 'small_task'
  if (params.publish) {
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.bgi'
  }

  input:
    tuple val(filename), file(bgen_file), file(bgi_file), val(chrom)

  output:
    tuple val(filename), path(bgen_file), path(bgi_file), file("${bgen_file.baseName}.sample"), val(chrom)

script:
"""
qctool -g ${bgen_file} -os ${bgen_file.baseName}.sample
"""
}