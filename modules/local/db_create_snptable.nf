process DB_CREATE_SNPTABLE {
    publishDir "${params.outdir}", mode: 'copy'

    label 'db_make'
    tag "${params.project}"

    input:
        tuple file(db_bcf_file), file(db_bcf_index)
        file(annotations)
        file(update_snps_sql)
        file(min_header)
        
    output:
        path "snps_index-${workflow.sessionId}.db", emit: snps_db

    script:
    def echtvar_sources = annotations.collect{ "-e $it" }.join(" ")
    """
    echo -e "\n== create snp vcf =="
    time (cat $min_header && echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" && bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t100\tPASS\t.\tGT\t.\n" $bcf_file) | bgzip -@${task.cpus} -c > snp.table.vcf.gz

    echo -e "\n== index snp table =="
    tabix -p vcf --csi snp.table.vcf.gz

    echo -e "\n== annotate snps =="
    echtvar anno $echtvar_sources snp.table.vcf.gz snp.table.anno.vcf.gz
    
    echo -e "\n== create snp table =="
    bcftools query -f "%CHROM\t%POS\t%CHROM-%POS-%REF-%ALT\t%REF\t%ALT\trs%dbsnp_rs\n" snp.table.anno.vcf.gz > snp.table

    echo -e "\n== store db =="
    cat $update_snps_sql | sqlite3 snps_index-${workflow.sessionId}.db
    """    
}