process CLUMP_RESULTS {
    lable 'plink_clump'
    publishDir "${params.outdir}/results/clump", mode: 'copy'

    input:
        tuple file(imputed_bgen)
        tuple val(phenotype), file(pheno_result_gz)

    output:
        tuple file("${prenotype}.clump"), file("${phenotype}.clump.range")

    script:
    """
    zcat $pheno_result_gz | awk '{OFS="\t"}; NR == 1 {print \$0, "PVAL"}; NR > 1 {print \$0, 10^(-\$13)}' > regenie.pval
    plink \
        --memory ${task.memory.toMega()} \
        --bgen $imputed_bgen \
        --clump regenie.pval \
        --clump-p1 1e-5 \
        --clump-kb 500 \
        --clump-r2 0.5 \
        --clump-snp-field ID \
        --clump-field PVAL \
        --clump-range glist-hg19 \
        --clump-range-border 20 \
        --out ${phenotype}
    """
}

