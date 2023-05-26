process MAKE_VARIANTS_CHUNKS {
    if (params.publish) {
        publishDir "${params.outdir}", mode: 'copy', pattern: "${filename}.GWAS-chunks.txt"
    }

    input:
        tuple val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom), file(snplist_file)

    output:
        tuple val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom), file("${filename}.GWAS-chunks.txt")

    script:
    def pos_idx = snplist_file.extension == "pvar" ? 2 : 4
    def chromosomes_list = chrom == "ONE_FILE" ? params.chromosomes.join(" ") : "$chrom"
    """
    grep -v "#"  $snplist_file | awk '{print \$1, \$${pos_idx} >> \$1".snps"}'
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
    if (params.publish) {
        publishDir "${params.outdir}", mode: 'copy', pattern: "${filename}.rarevar-chunks.chr*"
    }

    input:
        tuple val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom), file(set_list_file)
        //set list file according to regenie specs: gene_name, chrom, start_pos, list of vars

    output:
        tuple val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom), file("${filename}.rarevar-chunks.chr*")

    script:
    def chromosomes_list = chrom == "ONE_FILE" ? params.chromosomes.join(" ") : "$chrom"
    """
    grep -v "#"  $set_list_file | awk '{print \$1 >> \$2".genes"}'
    for c in $chromosomes_list
    do 
        if [ -f \$c.genes ]
        then
            split -l ${params.step2_rarevar_chunk_size} -d 1.genes ${filename}.rarevar-chunks.chr\${c}.
        fi
    done
    """
}