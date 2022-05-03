process DB_MAKE_VCF_PHENO {
    module 'htslib-tools/1.14'
    tag "$phenotype"

    input:
        tuple val(phenotype), file(regenie_results_gz)
        file(min_header)
        file(models_table)
        
    output:
        tuple file("${phenotype}.vcf.gz"), file("${phenotype}.vcf.gz.csi"), emit: pheno_db

    script:
    """
    cp $min_header header.txt
    echo "##MODEL=<ID=1,pheno=\\"$phenotype\\",vartype=\\"\$(grep -m1 $phenotype $models_table | cut -f3)\\",description=\\"\$(grep -m1 $phenotype $models_table | cut -f2)\\">" >> header.txt

    (cat header.txt && zcat $f \
    | awk -F"\t" 'BEGIN {OFS="\t"; print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "1"}; {OFS="\t"}; NR>1 && \$13 > ${params.pval_threshold} {if (\$13 > 7.301) GT="1/1"; else if (\$13 > 5) GT="0/1"; else GT="0/0"}; NR==1 {if (\$10 == "BETA") VARTYPE="Q"; else VARTYPE="B"}; NR > 1 {print \$1,\$2,\$3,\$4,\$5,"100","PASS",".","GT:INFO:LOGP:TYPE:EFFECT:SE:N:AF", GT":"\$7":"\$13":"VARTYPE":"\$10":"\$11":"\$8":"\$6}') \
    | bgzip -@${task.cpus} -c > ${phenotype}.vcf.gz
    
    tabix -p vcf --csi ${phenotype}.vcf.gz
    """    
}

process DB_MAKE_VCF_SINGLE {
    publishDir "${params.outdir}", mode: 'copy'
    
    module 'htslib-tools/1.14:bcftools/1.15'
    tag "$phenotype"

    input:
        tuple val(phenotype), file(regenie_results_gz)
        file(min_header)
        file(models_table)
        
    output:
        tuple file('gwas_db.bcf'), file('gwas_db.bcf.csi'), emit: pheno_db

    script:
    """
    cp $min_header header.txt
    echo "##MODEL=<ID=1,pheno=\\"$phenotype\\",vartype=\\"\$(grep -m1 $phenotype $models_table | cut -f3)\\",description=\\"\$(grep -m1 $phenotype $models_table | cut -f2)\\">" >> header.txt

    (cat header.txt && zcat $f \
    | awk -F"\t" 'BEGIN {OFS="\t"; print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "1"}; {OFS="\t"}; NR>1 && \$13 > ${params.pval_threshold} {if (\$13 > 7.301) GT="1/1"; else if (\$13 > 5) GT="0/1"; else GT="0/0"}; NR==1 {if (\$10 == "BETA") VARTYPE="Q"; else VARTYPE="B"}; NR > 1 {print \$1,\$2,\$3,\$4,\$5,"100","PASS",".","GT:INFO:LOGP:TYPE:EFFECT:SE:N:AF", GT":"\$7":"\$13":"VARTYPE":"\$10":"\$11":"\$8":"\$6}') \
    | bgzip -@${task.cpus} -c > ${phenotype}.vcf.gz
    
    bcftools view -Ob -o gwas_db.bcf ${phenotype}.vcf.gz
    bcftools index --csi gwas_db.bcf
    """    
}