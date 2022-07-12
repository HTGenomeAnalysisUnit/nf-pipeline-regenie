workflow CLUMP_RESULTS {
    take:
        results //tuple val(phenotype), file(pheno_results_gz)
        genes_interval_hg19 //file
        genes_interval_hg38 //file
        imputed_bgen //tuple(imputed_bgen_file.baseName, imputed_bgen_file, imputed_bgen_idx)
        sample_file //file

    main:
        if (params.ld_panel == 'NO_LD_FILE') {
            CONVERT_TO_BED(imputed_bgen, sample_file)
            chromosomes_ch = Channel.fromList(params.chromosomes)
            ld_panel_ch = chromosomes_ch.combine(CONVERT_TO_BED.out)
        } else {
            ld_panel_ch = Channel.fromFilePairs("${params.ld_panel.replace('{CHROM}','*')}.{bed,bim,fam}", size:3)
        }
        clump_input_ch = results.combine(ld_panel_ch)
        PLINK_CLUMPING(clump_input_ch, genes_interval_hg19, genes_interval_hg38)
        merge_input_ch = PLINK_CLUMPING.out.groupTuple()
        MERGE_CLUMP_RESULTS(merge_input_ch)
    
    emit:
        best_loci = MERGE_CLUMP_RESULTS.out.annotloci
}

process CONVERT_TO_BED {
    label 'process_plink2'
    
    input:
        tuple val(bgen_prefix), file(imputed_bgen_file), file(dummy_sample)
        path sample_file

    output:
        path "${bgen_prefix}.{bed,bim,fam}"

    script:
    def bgen_sample = sample_file.name != 'NO_SAMPLE_FILE' ? "--sample $sample_file" : ''
    """
    plink2 \
    --bgen imputed_bgen_file ref-first \
    --keep samples_ld.tsv \
    --make-bed \
    --memory ${task.memory.toMega()} \
    --threads ${task.cpus} \
    $bgen_sample \
    --out ${bgen_prefix} 
    """
}

process PLINK_CLUMPING {
    publishDir "${params.outdir}/logs/${phenotype}_clump", mode: 'copy', pattern: '*.log'
    label 'plink_clump'
    tag "${phenotype}_${chrom}"

    input:
        tuple val(phenotype), file(pheno_results_gz), val(chrom), file(bfile)
        path genes_hg19, stageAs: 'hg19_genes'
        path genes_hg38, stageAs: 'hg38_genes'

    output:
        tuple val(phenotype), file("${chrom}.clumped"), file("${chrom}.clumped.ranges")

    script:
    def bfile_prefix = bfile[0].simpleName
    def genes_ranges = params.genotypes_build == 'hg19' ? "hg19_genes" : "hg38_genes"
    """
    touch ${chrom}.clumped
    touch ${chrom}.clumped.ranges

    zcat $pheno_results_gz | awk '{OFS="\t"}; NR == 1 {print \$0, "PVAL", "SNPID"}; NR > 1 {print \$0, 10^(-\$13), \$3"_"\$4"_"\$5}' > regenie.pval
    plink \
        --memory ${task.memory.toMega()} \
        --bfile $bfile_prefix \
        --chr $chrom \
        --clump regenie.pval \
        --clump-p1 ${params.clump_p1} \
        --clump-p2 ${params.clump_p2} \
        --clump-kb ${params.clump_kb} \
        --clump-r2 0.5 \
        --clump-snp-field SNPID \
        --clump-field PVAL \
        --clump-range $genes_ranges \
        --clump-range-border ${params.annotation_interval_kb} \
        --out ${chrom}
    """
}

process MERGE_CLUMP_RESULTS {
    publishDir "${params.outdir}/results/toploci", mode: 'copy'
    label 'merge_clump'
    tag "${phenotype}"

    input:
        tuple val(phenotype), file(chromosome_clump), file(chromosome_ranges)

    output:
        path "${phenotype}.toploci.tsv", emit: toploci
        path "${phenotype}.toploci.annot.tsv", emit: annotloci

    script:
    """
    echo -e "CHR\tF\tSNP\tBP\tP\tTOTAL\tNSIG\tS05\tS01\tS001\tS0001\tSP2" > ${phenotype}.toploci.tsv
    for f in *.clumped
    do
        tail -n+2 \$f | tr -s " " "\\t" | sed 's/^\\t//g' >> toploci.tsv
    done
    sed '/^\$/d' toploci.tsv | sort -k5,5g >> ${phenotype}.toploci.tsv

    echo -e "CHR\tSNP\tP\tN\tPOS\tKB\tRANGES" > ${phenotype}.toploci.annot.tsv
    do
        tail -n+2 \$f | tr -s " " "\\t" | sed 's/^\\t//g' >> toploci.annot.tsv
    done
    sed '/^\$/d' toploci.annot.tsv | sort -k3,3g >> ${phenotype}.toploci.annot.tsv
    """
}