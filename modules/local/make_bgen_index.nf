/*
workflow CHECK_BGEN_INDEX {
  take:
    imputed_bgen //path(bgen), path(bgi), path(sample)

  main:
    //if(imputed_bgen[1].exists()) {
    //  output_channel = imputed_bgen
    //} else {
      MAKE_BGEN_INDEX(imputed_bgen)
      output_channel = MAKE_BGEN_INDEX.out
    //}

  emit:
    bgen_with_index = output_channel
}
*/

process CHECK_BGEN_INDEX {
  label 'process_bgenix'
  if (params.publish) {
    publishDir "${params.outdir}/bgen", mode: 'copy', pattern: '*.bgi'
  }

  input:
    tuple val(bgen_basename), path(imputed_bgen_file), path(bgi_index)

  output:
    tuple val(bgen_basename), path("${imputed_bgen_file}"), path("${imputed_bgen_file}.bgi")

script:
def has_bgi = bgi_index.exists() ? "1" : "0"
"""
if [ "$has_bgi" == "0" ]
then
  rm $bgi_index
  bgenix -g ${imputed_bgen_file} -index
fi
"""
}