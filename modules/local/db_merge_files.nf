process DB_MERGE {
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${params.project}_db.${filename.extension}"}

    tag "${params.project}"
    
    module 'bcftools/1.15'

    input:
        tuple file(bcf_files), file(bcf_indexes)
        file(min_header)
        
    output:
        tuple file("gwas_db-${workflow.sessionID}-${task.index}.bcf"), file("gwas_db-${workflow.sessionID}-${task.index}.bcf.csi"), emit: gwas_db

    script:
    """
    n=0
    for f in *.bcf
    do
        echo "== \$f =="
        echo "\$f" >> filelist
        bcftools view -h \$f | grep "##contig" >> contigs.txt

        while read sample
        do
            n=\$(( \$n + 1 ))
            bcftools view -h \$f | grep "##MODEL" | grep -w "ID=\$sample" | sed "s/ID=\$sample/ID=\$n/" >> new_models.txt
        done <<<\$(bcftools view -h \$f | grep -m1 "#CHROM" | cut -f10- | tr "\t" "\n")
    done

    sort -u contigs.txt > contigs.unique.txt
    cat $min_header new_models.txt contigs.unique.txt > header.txt

    echo "tot samples: \$n"
    for s in \$(seq 1 \$n)
    do
        samples="\$samples\t\$s"
    done
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\$samples" >> header.txt

    echo -e "\n== merge DBs =="
    bcftools merge --use-header header.txt --threads ${task.cpus} -m none --no-version -l filelist -Ob -o gwas_db.bcf

    echo -e "\n== index =="
    bcftools index --csi gwas_db-${workflow.sessionID}-${task.index}.bcf
    """    
}