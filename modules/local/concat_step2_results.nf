process CONCAT_STEP2_RESULTS {
    publishDir "${params.outdir}/results", mode: 'copy'
    
    input:
        tuple val(filename), file(regenie_gz)
    
    output:
        path "*regenie.gz", emit: regenie_results_gz
        path "*regenie.gz.tbi", emit: regenie_results_tbi
    
    script:
    """
    #make tmp dir for sort 
    mkdir tmp_sort

    #get list of phenos from regenie output file names
    for f in *.regenie.gz
    do 
        pheno=\${f%%.regenie.gz}
        pheno=\${pheno##*_}
        echo "\$pheno" >> phenos.tmp
    done
    sort -u phenos.tmp > phenos.list
    
    #concat chunk files for each pheno
    while read -r pheno
    do
        for f in *\${pheno}.regenie.gz;
        do
            zcat \$f | tail -n+2 >> ${filename}_\${pheno}.tmp
        done
        headerfile=\$(ls *\${pheno}.regenie.gz | head -1)
        zcat \$headerfile | head -1 | sed 's/ /\t/g' > header.txt
        (cat header.txt && sed 's/ /\t/g' ${filename}_\${pheno}.tmp | sort -k1,1V -k2,2n -T tmp_sort) | bgzip -c > \${pheno}.regenie.gz 
        tabix -f -b 2 -e 2 -s 1 -S 1 \${pheno}.regenie.gz
    done < phenos.list
    
    #remove chunk files
    rm *_${filename}_*.regenie.gz
    """
}