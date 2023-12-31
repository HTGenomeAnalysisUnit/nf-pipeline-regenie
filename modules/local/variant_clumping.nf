process CONVERT_TO_BED {
    label 'process_plink2'
    
    input:
        tuple val(chrom), file(bgen_pgen), file(bgi_pvar), file(sample_psam)

    output:
        tuple val(chrom), file("${fileprefix}.bed"), file("${fileprefix}.bim"), file("${fileprefix}.fam")

    script:
    def bgen_sample = params.genotypes_imputed_format == 'bgen' ? "--sample $sample_psam" : ''
    def format = params.genotypes_imputed_format == 'pgen' ? 'pfile' : 'bgen'
    def extension = params.genotypes_imputed_format in ['bgen','vcf','bcf'] ? '.bgen ref-first' : ''
    fileprefix = bgen_pgen.baseName
    """
    plink2 \
    --${format} ${fileprefix}${extension} \
    --make-bed \
    --memory ${task.memory.toMega()} \
    --threads ${task.cpus} \
    $bgen_sample \
    --out ${fileprefix}
    """
}

process MERGE_BED_DATASET {
    label 'process_plink2'
    
    input:
        file(bed_files)
        file(bim_files)
        file(fam_files)

    output:
        tuple file("ld_panel_merged.bed"), file("ld_panel_merged.bim"), file("ld_panel_merged.fam")

    script:
    def merge_list = bed_files.collect{ "${it.baseName}" }.join(',')
    """
    echo "${merge_list}" | tr "," "\\n" | sort -u > files_to_merge.list

    plink \
    --merge-list files_to_merge.list \
    --make-bed \
    --memory ${task.memory.toMega()} \
    --threads ${task.cpus} \
    --out ld_panel_merged
    """
}

process PLINK_CLUMPING {
    publishDir {"${params.logdir}/${project_id}/logs/${phenotype}_clump"}, mode: 'copy', pattern: '*.log'
    label 'plink_clump'
    tag "${project_id}_${phenotype}_${chrom}"

    input:
        tuple val(project_id), val(phenotype), file(pheno_results_gz), val(chrom), file(bed), file(bim), file(fam)
        path genes_hg19, stageAs: 'hg19_genes'
        path genes_hg38, stageAs: 'hg38_genes'

    output:
        tuple val(project_id), val(phenotype), path("${output_prefix}.clumped"), path("${output_prefix}.clumped.ranges"), emit: chrom_clump_results
        path "${output_prefix}.log", emit: logs

    script:
    def bfile_prefix = bed.baseName
    def genes_ranges = params.genotypes_build == 'hg19' ? "hg19_genes" : "hg38_genes"
    def target_chrom = chrom != 'NO_SPLIT' ? "--chr ${chrom}" : ''
    output_prefix = chrom != 'NO_SPLIT' ? "${chrom}" : "${bfile_prefix}"
    """
    touch ${output_prefix}.clumped
    touch ${output_prefix}.clumped.ranges

    zcat $pheno_results_gz | awk '{OFS="\t"}; NR == 1 {print \$0, "PVAL"}; NR > 1 {print \$0, 10^(-\$13)}' > regenie.pval
    plink \
        --memory ${task.memory.toMega()} \
        --bfile $bfile_prefix \
        $target_chrom \
        --clump regenie.pval \
        --clump-p1 ${params.clump_p1} \
        --clump-p2 ${params.clump_p2} \
        --clump-kb ${params.clump_kb} \
        --clump-r2 0.5 \
        --clump-snp-field ID \
        --clump-field PVAL \
        --clump-range $genes_ranges \
        --clump-range-border ${params.annotation_interval_kb} \
        --out ${output_prefix}
    """
}

process MERGE_CLUMP_RESULTS {
    publishDir {"${params.outdir}/${project_id}/results/gwas/toploci"}, mode: 'copy'
    label 'merge_clump'
    tag "${project_id}_${phenotype}"

    input:
        tuple val(project_id), val(phenotype), file(chromosome_clump), file(chromosome_ranges)

    output:
        tuple val(project_id), val(phenotype), path("${phenotype}.toploci.tsv"), emit: toploci
        tuple val(project_id), val(phenotype), path("${phenotype}.toploci.annot.tsv"), emit: annotloci

    script:
    """
    echo -e "CHR\tF\tSNP\tBP\tP\tTOTAL\tNSIG\tS05\tS01\tS001\tS0001\tSP2" > ${phenotype}.toploci.tsv
    for f in *.clumped
    do
        tail -n+2 \$f | tr -s " " "\\t" | sed 's/^\\t//g' >> toploci.tsv
    done
    sed '/^\$/d' toploci.tsv | sort -k5,5g >> ${phenotype}.toploci.tsv

    echo -e "CHR\tSNP\tP\tN\tPOS\tKB\tRANGES" > ${phenotype}.toploci.annot.tsv
    for f in *.clumped.ranges
    do
        tail -n+2 \$f | tr -s " " "\\t" | sed 's/^\\t//g' >> toploci.annot.tsv
    done
    sed '/^\$/d' toploci.annot.tsv | sort -k3,3g >> ${phenotype}.toploci.annot.tsv
    """
}
