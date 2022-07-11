process MAKE_SNPLIST {
    label 'process_bgenix'
    if (params.publish) {
        publishDir "${params.outdir}/snplist", mode: 'copy'
    }

    input:
        tuple val(filename), path(bgen_file), path(bgi_index)

    output:
        file "${filename}.snplist"

    script:
    """
    bgenix -g $bgen_file -list | tail -n+3 | cut -f3,4 > ${filename}.snplist
    """
}