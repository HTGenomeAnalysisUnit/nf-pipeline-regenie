//Given a bgen file, use bgenix to extract the snp list and save it as a tab delimited file
//The output file is a 6 column file using the same format as the .bim file
process MAKE_SNPLIST {
    label 'split_data'
    
    if (params.publish) {
        publishDir "${params.outdir}", mode: 'copy'
    }

    input:
        tuple val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom)

    output:
        tuple val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom), file("${filename}.snplist")

    script:
    """
    bgenix -g $bed_bgen_pgen -list | tail -n+3 | awk '{OFS="\t"}; {print \$3, \$2, 0, \$4, \$6, \$7}' | sed '\$d' > ${filename}.snplist
    """
}