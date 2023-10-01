include { MAKE_SNPLIST          } from '../modules/local/make_snplist'  addParams(outdir: "${params.outdir}/bgen_dataset", publish: params.save_snplist)
include { MAKE_VARIANTS_CHUNKS  } from '../modules/local/make_chunks'   addParams(outdir: "${params.outdir}/step2_chunks", chromosomes: params.chromosomes, publish: params.save_chunks_file, snplist_type: params.input_format)
include { MAKE_GENES_CHUNKS     } from '../modules/local/make_chunks'   addParams(outdir: "${params.outdir}/step2_chunks", chromosomes: params.chromosomes, publish: params.save_chunks_file)

workflow SPLIT_GWAS_DATA_WF {
    take:
    genotypes_files //[val(filename), file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom)]
                    // chrom = "ONE_FILE" if we don't split by chromosome

    main:
    if (params.input_format == 'bgen') {
    //==== MAKE A SNPLIST IF INPUT IS BGEN ====
        // If a .snplist file is provided we use it
        if (params.snplist_file) {
            snplist_file = file(params.snplist_file, checkIfExists: true)
            genotypes_files
                .map{ tuple(it[0], it[1], it[2], it[3], it[4], snplist_file) }
        } else if (params.snplist_folder) {
        // If a folder is provided we use the .snplist files inside
            genotypes_files
                .map{ tuple(it[0], it[1], it[2], it[3], it[4], file("${params.snplist_folder}/${it[1].baseName}.snplist")) }
                .branch {
                    found: it[5].exists()
                    missing: true
                }
                .set { snplist_files_ch }
        } else {
        //Search for .snplist file based on input bgen file
            genotypes_files
                .map{ tuple(it[0], it[1], it[2], it[3], it[4], file("${it[1].parent}/${it[1].baseName}.snplist")) }
                .branch {
                    found: it[5].exists()
                    missing: true
                }
            .set { snplist_files_ch }
        }

        //Make a snplist file if it doesn't exist
        new_snplist_ch = MAKE_SNPLIST(snplist_files_ch.missing.map { tuple(it[0], it[1], it[2], it[3], it[4]) })
        genotypes_snplist_ch = snplist_files_ch.found.mix(new_snplist_ch)
    } else {
    //==== WITH BED / PGEN USE BIM / PVAR DIRECTLY ====
        genotypes_snplist_ch = genotypes_files
            .map{ tuple(it[0], it[1], it[2], it[3], it[4], it[2]) }
    }

    //MAKE CHUNKS of N VARIANTS AS SPECIFIED IN PARAMS
    MAKE_VARIANTS_CHUNKS(genotypes_snplist_ch)
    genotypes_by_chunk_ch = MAKE_VARIANTS_CHUNKS.out
        .map{ tuple(it[0], it[1], it[2], it[3], it[4], it[5].readLines().flatten(), it[5].readLines().size()) }
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