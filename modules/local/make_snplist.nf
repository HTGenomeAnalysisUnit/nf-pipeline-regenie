process MAKE_SNPLIST {
    label 'process_plink2'
    if (params.publish) {
        publishDir "${params.outdir}/snplist", mode: 'copy'
    }

    input:
        tuple val(filename), path(bgen_file), path(bgi_index)

    output:
        file "${filename}.snplist"

    script:
    """
    plink2 \
        --bgen $bgen_file ref-first \
        --memory ${task.memory.toMega()} \
        --threads ${task.cpus} \
        --make-just-bim \
        --out temp
    
    cut -f1,4 temp.bim > ${filename}.snplist
    """
}