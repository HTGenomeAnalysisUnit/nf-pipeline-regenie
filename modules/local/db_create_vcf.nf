process DB_MAKE_VCF {
    input:
        tuple val(phenotype), file(regenie_results_gz)
        file(min_header)
        file(models_table)
        
    output:
        tuple file("${phenotype}.vcf.gz"), file("${phenotype}.vcf.gz.csi"), emit: pheno_vcf

    script:
    """
    f="$regenie_results_gz"
    
    cp $min_header header.txt
    echo "##MODEL=<ID=1,pheno=\\"$phenotype\\",vartype=\\"\$(grep -m1 $phenotype $models_table | cut -f3)\\",description=\\"\$(grep -m1 $phenotype $models_table | cut -f2)\\">" >> header.txt

    (cat header.txt && zcat $f \
    | awk -F"\t" 'BEGIN {OFS="\t"; print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "1"}; {OFS="\t"}; NR>1 && \$13 > ${params.pval_threshold} {if (\$13 > 7.301) GT="1/1"; else if (\$13 > 5) GT="0/1"; else GT="0/0"}; NR==1 {if (\$10 == "BETA") VARTYPE="Q"; else VARTYPE="B"}; NR > 1 {print \$1,\$2,\$3,\$4,\$5,"100","PASS",".","GT:INFO:LOGP:TYPE:EFFECT:SE:N:AF", GT":"\$7":"\$13":"VARTYPE":"\$10":"\$11":"\$8":"\$6}') \
    | bgzip -@4 -c > ${phenotype}.vcf.gz
    
    tabix -p vcf --csi ${phenotype}.vcf.gz
    """    
}