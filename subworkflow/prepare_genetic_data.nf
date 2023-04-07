workflow PREPARE_GENETIC_DATA {
    main:
    //==== PREPARE DATA FOR STEP2 ====
    if (params.genotypes_data_format == "vcf"){
        //convert vcf files to BGEN format and make BGI index
        genotypes_files =  Channel.fromPath("${params.genotypes_data}", checkIfExists: true)

        CONVERT_TO_BGEN (
            genotypes_files
        )
        genotypes_plink2_ch = CONVERT_TO_BGEN.out.genotypes_bgen
    } 
    if (params.genotypes_data_format == "bgen") {
        //Input is already BGEN 
        genotypes_bgen_and_index = Channel
            .fromPath(params.genotypes_data, checkIfExists: true)
            .map{ tuple(it.baseName, it, file("${it}.bgi")) }

        genotypes_bgen_and_sample = Channel
            .fromPath(params.genotypes_data, checkIfExists: true)
            .map{ tuple(it.baseName, file("${it.parent}/${it.baseName}.sample")) }

        //Check BGI index and make it if missing
        CHECK_BGEN_INDEX(genotypes_bgen_and_index)
        genotypes_plink2_ch = CHECK_BGEN_INDEX.out.join(genotypes_bgen_and_sample)
    }

    //Right now creating snplist from BGEN is time consuming for large data
    snplist_ch = Channel.empty()
    if (params.genotypes_data_format in ['bgen', 'vcf']) {
        if (params.step2_split_by == 'chunk') { //If we split by chunk we need snplist
            if (params.variants_snplist) { //If snplist is provided
                snplist_ch = Channel.fromPath("${params.genotypes_data}", checkIfExists: true)
                .map{ it -> return tuple(it.baseName, file("${it.parent}/${it.baseName}.snplist", checkIfExists: true)) }
            } else { //If snplist is not provided we make one from BGEN
                snplist_ch = MAKE_SNPLIST(genotypes_plink2_ch)
            }
        }
    } else {
        //If input is bed/bim/fam we can use the bim file as snplist
        genotypes_plink2_ch = Channel.fromFilePairs("${params.genotypes_data}.{bed,bim,fam}", size: 3)
        snplist_ch = tuple("${params.genotypes_data}", file("${params.genotypes_data}.snplist", checkIfExists: true))
    }

    emit:
        snplist_out = snplist_ch
        processed_genotypes_out = genotypes_plink2_ch
}