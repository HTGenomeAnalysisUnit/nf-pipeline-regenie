process CONCAT_STEP2_RESULTS {
    publishDir "${params.outdir}/${project_id}/results/gwas", mode: 'copy'
    
    input:
        tuple val(project_id), val(pheno), file(regenie_gz)
    
    output:
        tuple val(project_id), val(pheno), path("${pheno}.${suffix}.regenie.gz"), emit: regenie_results_gz
        path "${pheno}.${suffix}.regenie.gz.tbi", emit: regenie_results_tbi
    
    script:
    def n_head_lines = params.rarevar_results ? 2 : 1
    suffix = params.rarevar_results ? "rarevars" : "gwas"
    """
    mkdir -p tmp_sort

    for f in *${pheno}.regenie.gz;
    do
        zcat \$f | tail -n+${n_head_lines + 1} >> ${pheno}.tmp
    done
    
    headerfile=\$(ls *${pheno}.regenie.gz | head -1)

    zcat \$headerfile | head -n ${n_head_lines} | sed 's/ /\t/g' > header.txt
    (cat header.txt && sed 's/ /\t/g' ${pheno}.tmp | sort -S ${task.memory.toGiga() - 1}G -k1,1V -k2,2n -T tmp_sort) | bgzip -c > ${pheno}.${suffix}.regenie.gz
    
    # Check sort exit status, this is need to properly capture OOM errors on sort by returning the actual exit status
    # Without this check, the exit status of sort is always 0 due to pipe
    sort_exit=\${PIPESTATUS[0]}
    if [[ \$sort_exit -ne 0 ]]
    then 
        echo "ERROR CODE: \$sort_exit"
        exit \$sort_exit
    fi

    tabix -f -b 2 -e 2 -s 1 -S 1 ${pheno}.${suffix}.regenie.gz
    """
}
