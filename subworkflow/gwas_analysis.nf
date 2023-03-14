workflow GWAS_ANALYSIS {
    take:
        imputed_plink2_ch
        regenie_step1_out_ch
        phenotypes_file_validated
        covariates_file_validated
        regenie_log_parser_jar

    main:
    //==== REGENIE STEP 2 - COMMON VARS ====
    if (params.step2_split_by == 'chr') { 
        //PARALLELIZE BY CHROM
        chromosomes = Channel.of(params.chromosomes).flatten()
        processed_imputed_ch = imputed_plink2_ch.combine(chromosomes)
    } else if (params.step2_split_by == 'chunk') {
        //PARALLELIZE BY CHUNK
        MAKE_CHUNKS(snplist_ch)
        chunks_ch = MAKE_CHUNKS.out
        .map{ it -> return tuple(it[0], it[1].split('\n').flatten()) }
        .transpose(by: 1)
        processed_imputed_ch = imputed_plink2_ch.combine(chunks_ch, by: 0)
    } else {
        //WORK ON THE WHOLE FILE
        processed_imputed_ch = imputed_plink2_ch
        .map{ it -> return tuple(it[0], it[1], it[2], it[3], "SINGLEFILE")}
    }

    REGENIE_STEP2 (
        regenie_step1_out_ch.collect(),
        processed_imputed_ch,
        phenotypes_file_validated,
        covariates_file_validated
    )

    REGENIE_LOG_PARSER_STEP2 (
        REGENIE_STEP2.out.regenie_step2_out_log.first(),
        regenie_log_parser_jar
    )

    /*
    if (params.step2_split_by == 'chr' || params.step2_split_by == 'chunk') {
        //concat by chromosome results into a single result file per pheno
        concat_input_ch = REGENIE_STEP2.out.regenie_step2_out.transpose()
        .map{ it -> return tuple(it[2].baseName.replace("${it[0]}_${it[1]}_",'').replace('.regenie',''), it[2]) }
        .groupTuple().map{ it -> return tuple(it[0], it[1].flatten()) }
        
    } else {
        concat_input_ch = REGENIE_STEP2.out.regenie_step2_out
        .transpose()
        .map{ it -> return tuple(it[2].baseName.replace("${it[0]}_",'').replace('.regenie',''), it[1]) }
        .groupTuple().map{ it -> return tuple(it[0], it[1].flatten()) }
    }
    */

    concat_input_ch = REGENIE_STEP2.out.regenie_step2_out.transpose()
        .map{ it -> return tuple(it[2].baseName.replace("${it[0]}_${it[1]}_",'').replace('.regenie',''), it[2]) }
        .groupTuple().map{ it -> return tuple(it[0], it[1].flatten()) }
    CONCAT_STEP2_RESULTS(concat_input_ch)

    // generate a tuple of phenotype and corresponding result
    regenie_step2_by_phenotype = CONCAT_STEP2_RESULTS.out.regenie_results_gz

    //==== FILTER AND ANNOTATE TOP HITS ====
    FILTER_RESULTS (
        regenie_step2_by_phenotype
    )
    
    ANNOTATE_FILTERED (
        FILTER_RESULTS.out.results_filtered,
        genes_bed_hg19,
        genes_bed_hg38
    )
  
    //==== PERFORM VARIANT CLUMPING ====
    if (params.clumping) {
        if (params.ld_panel == 'NO_LD_FILE') {
            log.warn "No ld_panel provided, clumping will be performed using the whole genomic dataset"
        } 
        CLUMP_RESULTS(regenie_step2_by_phenotype, genes_ranges_hg19, genes_ranges_hg38, imputed_plink2_ch)
        clump_results_ch = CLUMP_RESULTS.out.best_loci
    } else {
        clump_results_ch = regenie_step2_by_phenotype.map { it -> return tuple(it[0], file('NO_CLUMP_FILE'))}
    }

    merged_results_and_annotated_filtered = regenie_step2_by_phenotype
        .join(ANNOTATE_FILTERED.out.annotated_ch, by: 0)
        .join(clump_results_ch, by: 0, remainder: true)

    //==== GENERATE HTML REPORTS ====
    if (params.make_report) {
        REPORT (
            merged_results_and_annotated_filtered,
            VALIDATE_PHENOTYPES.out.phenotypes_file_validated,
            gwas_report_template,
            VALIDATE_PHENOTYPES.out.phenotypes_file_validated_log,
            covariates_file_validated_log.collect(),
            regenie_step1_parsed_logs_ch.collect(),
            REGENIE_LOG_PARSER_STEP2.out.regenie_step2_parsed_logs
        )
    }

    emit:
        regenie_step2_out = regenie_step2_by_phenotype
        regenie_step2_filtered = FILTER_RESULTS.out.results_filtered
        regenie_step2_log = REGENIE_STEP2.out.regenie_step2_out_log.first()
}