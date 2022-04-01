process MAKE_CHUNKS {
    if (params.publish) {
        publishDir "${params.outdir}/step2_chunks", mode: 'copy', pattern: 'chunks.txt'
    }

    input:
        file(snps_list) //sorted tab-separated file with contig, snp pos
        val(chunk_size)

    output:
        path 'chunks.txt'

    script:
    """
    awk '{print \$0 >> \$1".snps"}' $snps_list
    for f in *.snps
    do 
        awk 'NR==1 || !(NR%10000); END { if (NR%10000) {print \$0} }' \$f \
        | paste - - \
        | awk 'NR==1 {print \$1":"\$2"-"\$4; start=\$4+1}; NR > 1 {print \$1":"start"-"\$2; print \$1":"\$2+1"-"\$4; start=\$4+1}' \
        > \$f.intervals
    done

    cat *.intervals > chunks.txt
    """
}