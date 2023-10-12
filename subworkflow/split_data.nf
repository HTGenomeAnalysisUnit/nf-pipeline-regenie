include { MAKE_VARIANTS_CHUNKS  } from '../modules/local/make_chunks'   addParams(outdir: "${params.outdir}/step2_chunks", chromosomes: params.chromosomes, publish: params.save_chunks_file, snplist_type: params.input_format)
include { MAKE_GENES_CHUNKS     } from '../modules/local/make_chunks'   addParams(outdir: "${params.outdir}/step2_chunks", chromosomes: params.chromosomes, publish: params.save_chunks_file)

workflow SPLIT_GWAS_DATA_WF {
    take:
    genotypes_files //[val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom)]
                    // chrom = "ONE_FILE" if we don't split by chromosome

    main:
    genotypes_snplist_ch = genotypes_files
        .map{ tuple(it[0], it[1], it[2], it[3], it[4], it[2]) }

    //MAKE CHUNKS of N VARIANTS AS SPECIFIED IN PARAMS
    MAKE_VARIANTS_CHUNKS(genotypes_snplist_ch, file("$projectDir/bin/get_intervals.sql", checkIfExists: true))
    count_chunks = MAKE_VARIANTS_CHUNKS.out
        .map{ it[5].readLines().size() } //count number of lines in each chunk
        .sum() //now sum total number of chunk for this analysis, this is used later for groupTuple
    
    genotypes_by_chunk_ch = MAKE_VARIANTS_CHUNKS.out
        .map{ tuple(it[0], it[1], it[2], it[3], it[4], it[5].readLines().flatten()) }
        .combine(count_chunks)
        .transpose(by: 5)

    emit:
        //val(filename), path(bed_bgen_pgen), path(bim_bgi_pvar), path(fam_sample_psam), val(chrom), path(chunk_coords), val(n_chunks)
        processed_genotypes = genotypes_by_chunk_ch
}

workflow SPLIT_RAREVARIANT_DATA_WF {
    take:
    genotypes_files //[val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom)]
                    // chrom = "ONE_FILE" if we don't split by chromosome
    
    main:
    genotypes_set_list = genotypes_files
        .map{ tuple(it[0], it[1], it[2], it[3], it[4], file(params.rarevar_set_list_file)) }
    MAKE_GENES_CHUNKS(genotypes_set_list)

    genotypes_by_chunk_ch = MAKE_GENES_CHUNKS.out
        .map{ tuple(it[0], it[1], it[2], it[3], it[4], it[5], it[5].size()) }
        .transpose(by: 5)

    emit:
        //val(filename), path(bed_bgen_pgen), path(bim_bgi_pvar), path(fam_sample_psam), val(chrom), path(chunk_genes), val(n_chunks)
        processed_genotypes = genotypes_by_chunk_ch
}