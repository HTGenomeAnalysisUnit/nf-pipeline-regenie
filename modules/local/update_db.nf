process UPDATE_DB {
    input:
        file(sql_lite_script)
        file(db_file)
        file(regenie_results)

    output:
        stdout

    script:
    """
    for f in *.regenie.gz
    do
        phenotype=\${f##.regenie.gz}
        tableid="${params.project}_\$phenotype"
        zcat \$f | awk -v pheno=\$phenotype '{OFS="\\t"}; \$13 >= 1.3 {print pheno, \$0}' > filtered.tsv
        
        cp $sql_lite_script \${tableid}.sql
        sed -i "s/@table_name@/\$tableid/g" \${tableid}.sql 
        sed -i "s/@filename@/filtered.tsv/g" \${tableid}.sql 
        cat \${tableid}.sql | sqlite3 $db_file
    done 
    echo "DONE"
    """
}