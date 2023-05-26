//Rmd template files
gwas_report_template = file("$projectDir/reports/gwas_report_template.Rmd",checkIfExists: true)
//rarevar_report_template = file("$projectDir/reports/rare_vars_report_template.Rmd",checkIfExists: true)

//Annotation files
if (params.genes_bed) {
  genes_bed_hg19 = file(params.genes_bed)
  genes_bed_hg38 = file(params.genes_bed)
} else {
  genes_bed_hg19 = file("$projectDir/genes/genes.GRCh37.1-23.sorted.bed", checkIfExists: true)
  genes_bed_hg38 = file("$projectDir/genes/genes.GRCh38.1-23.sorted.bed", checkIfExists: true)
}

if (params.genes_ranges) {
  genes_ranges_hg19 = file(params.genes_ranges)
  genes_ranges_hg38 = file(params.genes_ranges)
} else {
  genes_ranges_hg19 = file("$projectDir/genes/glist-hg19", checkIfExists: true)
  genes_ranges_hg38 = file("$projectDir/genes/glist-hg38", checkIfExists: true)
}

//Inclusion statements
include { FILTER_RESULTS              } from '../modules/local/filter_results'
include { ANNOTATE_FILTERED           } from '../modules/local/annotate_filtered'  addParams(outdir: "$outdir", annotation_interval_kb: params.annotation_interval_kb)
if (params.clumping) {
  include { CLUMP_RESULTS } from '../modules/local/clump_results' addParams(outdir: "$outdir")
}

workflow PROCESS_GWAS_RESULTS_WF {
	take:
    regenie_step2_by_phenotype
    imputed_plink2_ch

	main:
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
    //Ideally we can mix filtered_annotated channels for both rare vars and gwas and run reports for all
    html_reports_ch = Channel.empty()
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
        html_reports_ch = REPORT.out
    }

    emit:
    processed_results = merged_results_and_annotated_filtered
    html_reports = html_reports_ch
}

workflow PROCESS_RAREVAR_RESULTS_WF {
    take:
    regenie_step2_by_phenotype
    
	main:
	//==== FILTER AND ANNOTATE TOP HITS ====
    FILTER_RESULTS (
        regenie_step2_by_phenotype
    )
    
    /*
    //At the moment we don't provide any additional annotation for gene based results
    ANNOTATE_FILTERED (
        FILTER_RESULTS.out.results_filtered,
        genes_bed_hg19,
        genes_bed_hg38
    )
    */

    merged_results_and_annotated_filtered = regenie_step2_by_phenotype
        .join(FILTER_RESULTS.out.results_filtered, by: 0)
       //.map { tuple(it[0], it[1], it[2], "NO_CLUMP_FILE") }

    //==== GENERATE HTML REPORTS ====
    //Ideally we can mix filtered_annotated channels for both rare vars and gwas and run reports for all
    html_reports_ch = Channel.empty()
    if (params.make_report) {
        REPORT (
            merged_results_and_annotated_filtered,
            VALIDATE_PHENOTYPES.out.phenotypes_file_validated,
            rarevar_report_template,
            VALIDATE_PHENOTYPES.out.phenotypes_file_validated_log,
            covariates_file_validated_log.collect(),
            regenie_step1_parsed_logs_ch.collect(),
            REGENIE_LOG_PARSER_STEP2.out.regenie_step2_parsed_logs
        )
        html_reports_ch = REPORT.out
    }

    emit:
    processed_results = merged_results_and_annotated_filtered
    html_reports = html_reports_ch
}