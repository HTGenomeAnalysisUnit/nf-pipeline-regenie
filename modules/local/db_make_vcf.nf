process DB_MAKE_VCF_MULTI {
    label 'db_make'
    tag "$phenotype"

    input:
        tuple val(phenotype), file(regenie_results_gz), val(model), val(trait_type), val(genetic_model)
        file(min_header)
        
    output:
        tuple file("${phenotype}.vcf.gz"), file("${phenotype}.vcf.gz.csi"), emit: pheno_db

    script:
    """
    cp $min_header header.txt
    echo "##MODEL=<ID=1,pheno=\\"$phenotype\\",vartype=\\"$trait_type\\",genetic_model=\\"$genetic_model\\",description=\\"$model">" >> header.txt

    (cat header.txt && zcat $regenie_results_gz \
    | awk -F"\\t" 'BEGIN {OFS="\\t"; print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "1"}; {OFS="\\t"}; NR>1 && \$13 > ${params.pval_threshold} {if (\$13 > 7.301) GT="1/1"; else if (\$13 > 5) GT="0/1"; else GT="0/0"}; NR==1 {if (\$10 == "BETA") VARTYPE="Q"; else VARTYPE="B"}; NR > 1 {print \$1,\$2,\$3,\$4,\$5,"100","PASS",".","GT:INFO:LOGP:TYPE:EFFECT:SE:N:AF", GT":"\$7":"\$13":"VARTYPE":"\$10":"\$11":"\$8":"\$6}') \
    | bgzip -@${task.cpus} -c > ${phenotype}.vcf.gz
    
    tabix -p vcf --csi ${phenotype}.vcf.gz
    """    
}

process DB_MAKE_VCF_SINGLE {
    publishDir "${params.outdir}", mode: 'copy'
    
    label 'db_make'
    tag "$phenotype"

    input:
        tuple val(phenotype), file(regenie_results_gz), val(model), val(trait_type), val(genetic_model)
        file(min_header)
        
    output:
        tuple file("gwas_db-${workflow.sessionId}.bcf"), file("gwas_db-${workflow.sessionId}.bcf.csi"), emit: pheno_db

    script:
    """
    cp $min_header header.txt
    echo "##MODEL=<ID=1,pheno=\\"$phenotype\\",vartype=\\"$trait_type\\",genetic_model=\\"$genetic_model\\",description=\\"$model">" >> header.txt

    (cat header.txt && zcat $regenie_results_gz \
    | awk -F"\\t" 'BEGIN {OFS="\\t"; print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","1"}; {OFS="\\t"}; NR>1 && \$13 > ${params.pval_threshold} {if (\$13 > 7.301) GT="1/1"; else if (\$13 > 5) GT="0/1"; else GT="0/0"}; NR==1 {if (\$10 == "BETA") VARTYPE="Q"; else VARTYPE="B"}; NR > 1 {print \$1,\$2,\$3,\$4,\$5,"100","PASS",".","GT:INFO:LOGP:TYPE:EFFECT:SE:N:AF", GT":"\$7":"\$13":"VARTYPE":"\$10":"\$11":"\$8":"\$6}') \
    | bgzip -@${task.cpus} -c > ${phenotype}.vcf.gz
    
    bcftools view -Ob -o gwas_db-${workflow.sessionId}.bcf ${phenotype}.vcf.gz
    bcftools index --csi gwas_db-${workflow.sessionId}.bcf
    """    
}