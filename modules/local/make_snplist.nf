//Given a bgen file, use bgenix to extract the snp list and save it as a tab delimited file
//The output file is a 6 column file using the same format as the .bim file
process MAKE_SNPLIST {
    label 'process_bgenix'
    if (params.publish) {
        publishDir "${params.outdir}/snplist", mode: 'copy'
    }

    input:
        tuple val(filename), path(bgen_file), path(bgi_index), path(sample_file)

    output:
        tuple val(filename), file("${filename}.snplist")

    script:
    """
    bgenix -g $bgen_file -list | tail -n+3 | awk '{OFS="\t"}; {print \$3, \$2, 0, \$4, \$6, \$7}' | sed '\$d' > ${filename}.snplist
    """
}