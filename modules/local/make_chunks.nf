process MAKE_VARIANTS_CHUNKS {
    label 'small_task'
    
    if (params.publish) {
        publishDir "${params.outdir}", mode: 'copy', pattern: "*.GWAS-chunks.txt"
    }

    input:
        tuple val(filename), path(bed_bgen_pgen), path(bim_bgi_pvar), path(fam_sample_psam), val(chrom), path(snplist_file, stageAs: 'snplist.tsv')

    output:
        tuple val(filename), path(bed_bgen_pgen), path(bim_bgi_pvar), path(fam_sample_psam), val(chrom), path("${filename}.GWAS-chunks.txt")

    script:
    def pos_idx = params.snplist_type == "pgen" ? 2 : 4
    def chromosomes_list = chrom == "ONE_FILE" ? params.chromosomes.join(" ") : "$chrom"
    """
    grep -v "#"  snplist.tsv | awk '{print \$1, \$${pos_idx} >> \$1".snps"}'
    for c in $chromosomes_list
    do 
        if [ -f \$c.snps ]
        then
            awk 'NR == 1 {start=\$2}; NR > 1 && !(NR%${params.step2_gwas_chunk_size}) {print \$1":"start"-"\$2; start = \$2+1}; END { if (NR%${params.step2_gwas_chunk_size}) {print \$1":"start"-"\$2} }' \$c.snps > \$c.intervals 
        fi
    done

    cat *.intervals | sort -V > ${filename}.GWAS-chunks.txt
    """
}

process MAKE_GENES_CHUNKS {
    label 'small_task'
    
    if (params.publish) {
        publishDir "${params.outdir}", mode: 'copy', pattern: "*.rarevar-chunks.chr*"
    }

    input:
        tuple val(filename), path(bed_bgen_pgen), path(bim_bgi_pvar), path(fam_sample_psam), val(chrom), path(set_list_file)
        //set list file according to regenie specs: gene_name, chrom, start_pos, list of vars

    output:
        tuple val(filename), path(bed_bgen_pgen), path(bim_bgi_pvar), path(fam_sample_psam), val(chrom), path("${filename}.rarevar-chunks.chr*")

    script:
    def chromosomes_list = chrom == "ONE_FILE" ? params.chromosomes.join(" ") : "$chrom"
    """
    grep -v "#" $set_list_file | awk '{print \$1 >> \$2".genes"}'
    for c in $chromosomes_list
    do 
        if [ -f \$c.genes ]
        then
            split -l ${params.step2_rarevar_chunk_size} -d \${c}.genes ${filename}.rarevar-chunks.chr\${c}.
        fi
    done
    """
}