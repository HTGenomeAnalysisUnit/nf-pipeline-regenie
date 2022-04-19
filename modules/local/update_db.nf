process UPDATE_DB {
    input:
        file(save_results_sql)
        file(update_projects_sql)
        file(db_file)
        file(regenie_results)

    output:
        stdout

    script:
    def trait_type = params.phenotypes_binary_trait ? "binary" : "quantitative"
    """
    cp $update_projects_sql insert_new_projects.sql

    #store results
    for f in *.regenie.gz
    do
        phenotype=\${f%%.regenie.gz}
        tableid=\$(echo "${params.project}_\$phenotype" | tr "-" "_" | tr " " "_" | tr '.' '_')
        zcat \$f | awk -v pheno=\$phenotype -v runid="${params.project}" '{OFS="\\t"}; NR == 1 {print "RUN_ID","PHENO", \$0}; NR > 1 && \$13 >= ${params.db_pval_threshold} {print runid, pheno, \$0}' > filtered.tsv
        
        echo '("${params.project}", "'\$phenotype'", "${trait_type}", "${params.covariates_columns}", 500000, "${workflow.manifest.version}", "'\$USER'"),' >> insert_new_projects.sql

        cp $save_results_sql \${tableid}.sql
        sed -i "s/@table_name@/\$tableid/g" \${tableid}.sql 
        sed -i "s/@filename@/filtered.tsv/g" \${tableid}.sql 
        cat \${tableid}.sql | sqlite3 $db_file
    done

    #update projects table
    sed -i '$s/,/;/' insert_new_projects.sql
    cat insert_new_projects.sql | sqlite3 $db_file
   
    echo "DONE"
    """
}