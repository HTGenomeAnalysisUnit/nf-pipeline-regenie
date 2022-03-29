process CONCAT_STEP2_RESULTS {
    input:
        tuple val(filename), path(regenie_gz)
    
    output:
        tuple val(filename), path("${filename}*regenie.gz"), emit: regenie_step2_out
    
    script:
    """
    for f in *.regenie.gz; do pheno = \${f%%.regenie.gz}; pheno = \${pheno##_*} > phenos.list
    cat phenos.list | while read -r pheno
    do
        for f in ls *\${pheno}.regenie.gz;
        do
            zcat \$f | tail -n+2 >> ${filename}_\${pheno}.regenie
        done
        gzip ${filename}_\${pheno}.regenie
    done
    """
}