process CONCAT_STEP2_RESULTS {
    publishDir "${params.outdir}/results", mode: 'copy'
    
    input:
        tuple val(pheno), file(regenie_gz)
    
    output:
        tuple val(pheno), path("${pheno}.regenie.gz"), emit: regenie_results_gz
        path "${pheno}.regenie.gz.tbi", emit: regenie_results_tbi
    
    script:
    """
    mkdir tmp_sort

    for f in *${pheno}.regenie.gz;
    do
        zcat \$f | tail -n+2 >> ${pheno}.tmp
    done
    
    headerfile=\$(ls *${pheno}.regenie.gz | head -1)
    zcat \$headerfile | head -1 | sed 's/ /\t/g' > header.txt
    (cat header.txt && sed 's/ /\t/g' ${pheno}.tmp | sort -S ${task.memory.toGiga()}G -k1,1V -k2,2n -T tmp_sort) | bgzip -c > ${pheno}.regenie.gz 
    tabix -f -b 2 -e 2 -s 1 -S 1 ${pheno}.regenie.gz
    """
}