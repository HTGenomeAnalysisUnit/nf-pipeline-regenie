//Annotation files
if (params.genes_bed) {
  genes_bed_hg19 = file(params.genes_bed)
  genes_bed_hg38 = file(params.genes_bed)
} else {
  genes_bed_hg19 = file("$projectDir/genes/genes.GRCh37.1-23.${params.genes_group}.bed", checkIfExists: true)
  genes_bed_hg38 = file("$projectDir/genes/genes.GRCh38.1-23.${params.genes_group}.bed", checkIfExists: true)
}

if (params.genes_ranges) {
  genes_ranges_hg19 = file(params.genes_ranges)
  genes_ranges_hg38 = file(params.genes_ranges)
} else {
  genes_ranges_hg19 = file("$projectDir/genes/glist-hg19-${params.genes_group}", checkIfExists: true)
  genes_ranges_hg38 = file("$projectDir/genes/glist-hg38-${params.genes_group}", checkIfExists: true)
}

//Inclusion statements
include { FILTER_RESULTS    } from '../modules/local/filter_results'   addParams(outdir: "${params.outdir}", publish: params.publish_filtered)
include { ANNOTATE_FILTERED } from '../modules/local/annotate_filtered'  addParams(outdir: "${params.outdir}", annotation_interval_kb: params.annotation_interval_kb)
if (params.clumping) {
  include { CLUMP_RESULTS } from './clump_results' addParams(outdir: "${params.outdir}", logdir: "${params.logdir}", chromosomes: params.chromosomes)
}

workflow PROCESS_GWAS_RESULTS_WF {
	take:
    regenie_step2_by_phenotype //[project_id, phenotype, results_gz_file]
    processed_gwas_genotypes //[file_prefix, file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom)]

	main:
	//==== FILTER AND ANNOTATE TOP HITS ====
    FILTER_RESULTS ( regenie_step2_by_phenotype )
    
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
        CLUMP_RESULTS(regenie_step2_by_phenotype, genes_ranges_hg19, genes_ranges_hg38, processed_gwas_genotypes)
        clump_results_ch = CLUMP_RESULTS.out.best_loci
    } else {
        clump_results_ch = regenie_step2_by_phenotype.map { it -> return tuple(it[0], it[1], file('NO_CLUMP_FILE'))}
    }

    merged_results_and_annotated_filtered = regenie_step2_by_phenotype
        .join(ANNOTATE_FILTERED.out.annotated_ch, by: [0,1])
        .join(clump_results_ch, by: [0,1], remainder: true)

    emit:
    processed_results = merged_results_and_annotated_filtered
    //html_reports = html_reports_ch
}

workflow PROCESS_RAREVAR_RESULTS_WF {
    take:
    regenie_step2_by_phenotype //[project_id, phenotype, results_gz_file]
    
	main:
	//==== FILTER AND ANNOTATE TOP HITS ====
    FILTER_RESULTS ( regenie_step2_by_phenotype )
    
    /*
    //At the moment we don't provide any additional annotation for gene based results
    ANNOTATE_FILTERED (
        FILTER_RESULTS.out.results_filtered,
        genes_bed_hg19,
        genes_bed_hg38
    )
    */

    merged_results_and_annotated_filtered = regenie_step2_by_phenotype
        .join(FILTER_RESULTS.out.results_filtered, by: [0,1])
       //.map { tuple(it[0], it[1], it[2], "NO_CLUMP_FILE") }

    emit:
    processed_results = merged_results_and_annotated_filtered
}