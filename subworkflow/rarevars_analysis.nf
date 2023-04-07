workflow RARE_VARIANTS_ANALYSIS {
    take:
        imputed_plink2_ch
        regenie_step1_out_ch
        phenotypes_file_validated
        covariates_file_validated
        regenie_log_parser_jar

    main:

    //Read the gene list from the set list file and eventually split by gene
    if(params.step2_split_by == 'gene') {
        gene_list_ch = Channel.fromPath(params.rarevars_set_list_file)
            .splitCsv(header: false, sep: "\t").map{tuple(it[0], it[1])}

        // When input data has multiple files one per chromosome, we need to find the file corresponding to the gene
        // Here we expect the file name to contain a {CHROM} placeholder that will be replaced by the chromosome number		
        // SCENARIO 1: input data is a single file
        // SCENARIO 2: input data is a set of files, one per chromosome: use the {CHROM} placeholder in the file name
        // SCENARIO 3: input data is a generic set of files 

        // SCENARIO 1 accept split by chrom or gene
        // SCENARIO 2 only accept split by gene
        // SCENARIO 3 does not accept any split

        // tuple(filename.base, chrom(or all), bed_file, bim_file, fam_file), val(chunk/gene or all)
        hromosomes = params.chromosomes + params.sex_chromosomes.split(',')
    
        if (params.ld_panel =~ /\{CHROM\}/) 
            ld_panel_ch = Channel.of("${params.ld_panel}").combine(chromosomes.flatten())
            .map{ it[0].replace('{CHROM}', "${it[1]}") }
            .map{ tuple(file("${it}.bed", checkIfExist: true), file("${it}.bim", checkIfExist: true), file("${it}.fam", checkIfExist: true)) }
        else {
            ld_panel_ch = Channel.fromFilePairs("${params.ld_panel}.{bed,bim,fam}", size:3)
        }

        gene_list_ch = gene_list_ch.map{gene -> [gene, imputed_plink2_ch.find{it =~ /chr${gene}.*/}]}
    }

    //==== STEP2 FOR RARE VARIANTS ANALYSIS ====
    //we can use --extract-setlist to extract a specific gene in the set list file

    //==== CONCAT RESULTS ====

    //==== FILTER RESULTS AND ANNOTATE ====

    //==== GENERATE HTML REPORT ====

    emit:
}