process CONCAT_STEP2_RESULTS {
    input:
        tuple val(filename), file(regenie_gz)
    
    output:
        tuple val(filename), file("${filename}_*regenie.gz"), emit: regenie_step2_out
    
    script:
    """
    for f in *.regenie.gz; do pheno=\${f%%.regenie.gz}; pheno=\${pheno##*_}; echo "\$pheno" >> phenos.tmp; done
    sort -u phenos.tmp > phenos.list
    while read -r pheno
    do
        for f in *\${pheno}.regenie.gz;
        do
            zcat \$f | tail -n+2 >> ${filename}_\${pheno}.tmp
        done
        headerfile=\$(ls *\${pheno}.regenie.gz | head -1)
        zcat \$headerfile | head -1 > header.txt
        cat header.txt ${filename}_\${pheno}.tmp | gzip -c > ${filename}_\${pheno}.regenie.gz 
    done < phenos.list
    """
}