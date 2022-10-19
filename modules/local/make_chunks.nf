process MAKE_CHUNKS {
    if (params.publish) {
        publishDir "${params.outdir}/step2_chunks", mode: 'copy', pattern: 'chunks.txt'
    }

    input:
        tuple val(filename), file(snps_list) //sorted tab-separated file with contig, snp pos

    output:
        tuple val(filename), stdout

    script:
    def chromosomes_list = params.chromosomes.join(" ")
    """
    awk '{print \$0 >> \$1".snps"}' $snps_list
    for c in $chromosomes_list
    do 
        if [ -f \$c.snps ]
        then
            awk 'NR == 1 {start=\$2}; NR > 1 && !(NR%${params.step2_chunk_size}) {print \$1":"start"-"\$2; start = \$2+1}; END { if (NR%${params.step2_chunk_size}) {print \$1":"start"-"\$2} }' \$c.snps > \$c.intervals 
        fi
    done

    cat *.intervals | sort -V
    """
}