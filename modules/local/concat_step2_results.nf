process CONCAT_STEP2_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
        tuple val(pheno), file(regenie_gz)
    
    output:
        tuple val(pheno), path("${pheno}.${suffix}.regenie.gz"), emit: regenie_results_gz
        path "${pheno}.${suffix}.regenie.gz.tbi", emit: regenie_results_tbi
    
    script:
    def n_head_lines = params.rarevar_results ? 3 : 2
    suffix = params.rarevar_results ? "rarevars" : "gwas"
    """
    mkdir tmp_sort

    for f in *${pheno}.regenie.gz;
    do
        zcat \$f | tail -n+${n_head_lines} >> ${pheno}.tmp
    done
    
    headerfile=\$(ls *${pheno}.regenie.gz | head -1)
    zcat \$headerfile | head -1 | sed 's/ /\t/g' > header.txt
    (cat header.txt && sed 's/ /\t/g' ${pheno}.tmp | sort -k1,1V -k2,2n -T tmp_sort) | bgzip -c > ${pheno}.${suffix}.regenie.gz 
    tabix -f -b 2 -e 2 -s 1 -S 1 ${pheno}.${suffix}.regenie.gz
    """
}
