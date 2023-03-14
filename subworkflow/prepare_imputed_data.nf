workflow PREPARE_IMPUTED_DATA {
	take:

    main:
    //==== PREPARE IMPUTED DATA FOR STEP2 ====
    if (params.genotypes_imputed_format == "vcf"){
        //convert vcf files to BGEN format and make BGI index
        imputed_files =  Channel.fromPath("${params.genotypes_imputed}", checkIfExists: true)

        IMPUTED_TO_BGEN (
            imputed_files
        )
        imputed_plink2_ch = IMPUTED_TO_BGEN.out.imputed_bgen
    } else {
        //Input is already BGEN  
        imputed_bgen_and_index = Channel
            .fromPath(params.genotypes_imputed, checkIfExists: true)
            .map{ tuple(it.baseName, it, file("${it}.bgi")) }

        imputed_bgen_and_sample = Channel
            .fromPath(params.genotypes_imputed, checkIfExists: true)
            .map{ tuple(it.baseName, file("${it.parent}/${it.baseName}.sample")) }

        //Check BGI index and make it if missing
        CHECK_BGEN_INDEX(imputed_bgen_and_index)
        imputed_plink2_ch = CHECK_BGEN_INDEX.out.join(imputed_bgen_and_sample)
    }

    //If a snplist is not provided create one from BGEN
    //Right now this is time consuming for large data since need to convert format
    if (params.step2_split_by == 'chunk') {
        if (params.imputed_snplist && params.genotypes_imputed_format == "bgen") {
            snplist_ch = Channel.fromPath("${params.genotypes_imputed}", checkIfExists: true)
            .map{ it -> return tuple(it.baseName, file("${it.parent}/${it.baseName}.snplist", checkIfExists: true)) }
        } else {
            snplist_ch = MAKE_SNPLIST(imputed_plink2_ch)
        }
    }

    emit:
}