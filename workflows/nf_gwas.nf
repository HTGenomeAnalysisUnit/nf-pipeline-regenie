
requiredParams = [
    'project', 'genotypes_array',
    'genotypes_imputed', 'genotypes_build',
    'genotypes_imputed_format', 'phenotypes_filename',
    'phenotypes_columns', 'phenotypes_binary_trait',
    'regenie_test', 'step1_n_chunks', 'step2_split_by',
    'chromosomes'
]

for (param in requiredParams) {
    if (params[param] == null) {
      exit 1, "Parameter ${param} is required."
    }
}

if(params.outdir == null) {
  outdir = "${params.project}"
} else {
  outdir = "${params.outdir}/${params.project}"
}

if (params.master_log_dir != null) {
  master_log_dir = "${params.master_log_dir}"
} else {
  master_log_dir = "${outdir}"
}

project_id = params.project
phenotypes_array = params.phenotypes_columns.trim().split(',')

covariates_array= []
if(!params.covariates_columns.isEmpty()){
  covariates_array = params.covariates_columns.trim().split(',')
}

gwas_report_template = file("$projectDir/reports/gwas_report_template.Rmd",checkIfExists: true)

//Check scripts and resources
regenie_log_parser_java  = file("$projectDir/bin/RegenieLogParser.java", checkIfExists: true)
regenie_filter_java = file("$projectDir/bin/RegenieFilter.java", checkIfExists: true)
regenie_validate_input_java = file("$projectDir/bin/RegenieValidateInput.java", checkIfExists: true)

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

//Phenotypes
phenotypes_file = file(params.phenotypes_filename, checkIfExists: true)
phenotypes = Channel.from(phenotypes_array)

//Covariates
if (params.covariates_filename == 'NO_COV_FILE') {
  covar_tmp_file = file("${workflow.workDir}/NO_COV_FILE")
  covar_tmp_file.append('')
  covariates_file = file("$covar_tmp_file")
} else {
  covariates_file = file(params.covariates_filename)
}
if (params.covariates_filename != 'NO_COV_FILE' && !covariates_file.exists()){
  exit 1, "Covariate file ${params.covariates_filename} not found."
}


//Optional sample file
sample_file = file(params.regenie_sample_file)
if (params.regenie_sample_file != 'NO_SAMPLE_FILE' && !sample_file.exists()){
  exit 1, "Sample file ${params.regenie_sample_file} not found."
}

//Check specified test
if (params.regenie_test != 'additive' && params.regenie_test != 'recessive' && params.regenie_test != 'dominant'){
  exit 1, "Test ${params.regenie_test} not supported."
}

//Check imputed file format
if (params.genotypes_imputed_format != 'vcf' && params.genotypes_imputed_format != 'bgen'){
  exit 1, "File format ${params.genotypes_imputed_format} not supported."
}

//Array genotypes
Channel.fromFilePairs("${params.genotypes_array}.{bim,bed,fam}", size: 3).set {genotyped_plink_ch}

//Include modules
include { CACHE_JBANG_SCRIPTS         } from '../modules/local/cache_jbang_scripts'
include { VALIDATE_PHENOTYPES         } from '../modules/local/validate_phenotypes' addParams(outdir: "$outdir")
include { VALIDATE_COVARIATS          } from '../modules/local/validate_covariates' addParams(outdir: "$outdir")
include { IMPUTED_TO_BGEN             } from '../modules/local/imputed_to_plink2' addParams(outdir: "$outdir")
include { MAKE_BGEN_INDEX             } from '../modules/local/make_bgen_index' addParams(outdir: "$outdir", publish: params.save_bgen_index)
include { MAKE_SNPLIST                } from '../modules/local/make_snplist' addParams(outdir: "$outdir", publish: params.save_snplist)
include { PRUNE_GENOTYPED             } from '../modules/local/prune_genotyped' addParams(outdir: "$outdir")
include { QC_FILTER_GENOTYPED         } from '../modules/local/qc_filter_genotyped' addParams(outdir: "$outdir")
include { REGENIE_STEP1_SPLIT as REGENIE_STEP1 } from '../modules/local/regenie_step1' addParams(outdir: "$outdir", save_step1_predictions: params.save_step1_predictions, use_loocv: params.step1_use_loocv, niter: params.step1_niter)
include { REGENIE_LOG_PARSER_STEP1    } from '../modules/local/regenie_log_parser_step1'  addParams(outdir: "$outdir")
include { REGENIE_LOG_PARSER_STEP2    } from '../modules/local/regenie_log_parser_step2'  addParams(outdir: "$outdir")
include { FILTER_RESULTS              } from '../modules/local/filter_results'
include { MERGE_RESULTS_FILTERED      } from '../modules/local/merge_results_filtered'  addParams(outdir: "$outdir")
include { MERGE_RESULTS               } from '../modules/local/merge_results'  addParams(outdir: "$outdir")
include { ANNOTATE_FILTERED           } from '../modules/local/annotate_filtered'  addParams(outdir: "$outdir", annotation_interval_kb: params.annotation_interval_kb)
include { REPORT                      } from '../modules/local/report'  addParams(outdir: "$outdir")
include { CONCAT_STEP2_RESULTS        } from '../modules/local/concat_step2_results' addParams(outdir: "$outdir")
if (params.clumping) {
  include { CLUMP_RESULTS } from '../modules/local/clump_results' addParams(outdir: "$outdir")
}

if (params.step2_split_by == 'chunk') {
  include { MAKE_CHUNKS               } from '../modules/local/make_chunks.nf' addParams(publish: true, outdir: "$outdir", step2_chunk_size: params.step2_chunk_size, chromosomes: params.chromosomes)
  include { REGENIE_STEP2_BYCHUNK as REGENIE_STEP2     } from '../modules/local/regenie_step2' addParams(outdir: "$outdir", save_step2_logs: params.save_step2_logs)
} else if (params.step2_split_by == 'chr') {
  include { REGENIE_STEP2_BYCHR as REGENIE_STEP2     } from '../modules/local/regenie_step2' addParams(outdir: "$outdir", save_step2_logs: params.save_step2_logs)
}

//==== WORKFLOW ====
workflow NF_GWAS {   

    log_params = [ 
    'genotypes_array','genotypes_imputed','genotypes_imputed_format',
    'genotypes_build','imputed_snplist',
    'phenotypes_filename','phenotypes_columns','phenotypes_binary_trait',
    'covariates_filename','covariates_columns',
    'chromosomes',
    'regenie_test','annotation_min_log10p',
    'save_step1_predictions','outdir'
    ]
    log_params_string = []
    for (p in log_params) {
       log_params_string.add("$p : " + params[p])
    }

log.info"""\
==========================================================
  REGENIE GWAS - PROJECT ${params.project} - NF PIPELINE    
==========================================================
PROJECT ID      : ${params.project}

${log_params_string.join('\n')}
==========================================================
Please report issues to:
https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie
or contact: edoardo.giacopuzzi@fht.org
"""

    CACHE_JBANG_SCRIPTS (
        regenie_log_parser_java,
        regenie_filter_java,
        regenie_validate_input_java
    )

    VALIDATE_PHENOTYPES (
        phenotypes_file,
        CACHE_JBANG_SCRIPTS.out.regenie_validate_input_jar
    )

    if(covariates_file.exists() && covariates_file.name != 'NO_COV_FILE') {
        VALIDATE_COVARIATS (
          covariates_file,
          CACHE_JBANG_SCRIPTS.out.regenie_validate_input_jar
        )

        covariates_file_validated = VALIDATE_COVARIATS.out.covariates_file_validated
        covariates_file_validated_log = VALIDATE_COVARIATS.out.covariates_file_validated_log

   } else {

     // set covariates_file to default value
     covariates_file_validated = covariates_file
     covariates_file_validated_log = Channel.fromPath("NO_COV_LOG")

   }

    //convert vcf files to BGEN format and make BGI index
    if (params.genotypes_imputed_format == "vcf"){
      imputed_files =  channel.fromPath("${params.genotypes_imputed}", checkIfExists: true)

      IMPUTED_TO_BGEN (
          imputed_files
      )
      imputed_plink2_ch = IMPUTED_TO_BGEN.out.imputed_bgen
  
    } else {
    //Input is already BGEN  
      imputed_bgen_file = file(params.genotypes_imputed, checkIfExists: true)
      
      bgen_index_file = file("${params.genotypes_imputed}.bgi")
      //Make BGI index if missing otherwise set channel directly
      if (bgen_index_file.exists()) {
        channel.fromPath("${params.genotypes_imputed}")
        .map { tuple(it.baseName, it, file(it+".bgi")) }
        .set {imputed_plink2_ch}
      } else {
        MAKE_BGEN_INDEX(imputed_bgen_file)
        MAKE_BGEN_INDEX.out
        .map { tuple(it[0].baseName, it[0], it[1]) }
        .set {imputed_plink2_ch}
      }
    }

    //If a snplist is not provided create one from BGEN
    //Right now this is time consuming for large data since need to convert format
    if (params.step2_split_by == 'chunk') {
      if (params.imputed_snplist) {
        snplist_ch = file(params.imputed_snplist, checkIfExists: true)
      } else {
        snplist_ch = MAKE_SNPLIST(imputed_plink2_ch)
      }
    }

    QC_FILTER_GENOTYPED (
      genotyped_plink_ch, phenotypes_file
    )

    if(params.prune_enabled) {

      PRUNE_GENOTYPED (
          QC_FILTER_GENOTYPED.out.genotyped_filtered_files_ch
      )

      genotyped_final_ch = PRUNE_GENOTYPED.out.genotypes_pruned_ch

    } else {
      //no pruning applied, set QCed directly to genotyped_final_ch
      genotyped_final_ch = QC_FILTER_GENOTYPED.out.genotyped_filtered_files_ch
    }

    //==== REGENIE STEP 1 ====
    if (!params.regenie_premade_predictions){

      REGENIE_STEP1 (
          genotyped_final_ch,
          QC_FILTER_GENOTYPED.out.genotyped_filtered_snplist_ch,
          QC_FILTER_GENOTYPED.out.genotyped_filtered_id_ch,
          VALIDATE_PHENOTYPES.out.phenotypes_file_validated,
          covariates_file_validated
      )

      REGENIE_LOG_PARSER_STEP1 (
          REGENIE_STEP1.out.regenie_step1_out_log,
          CACHE_JBANG_SCRIPTS.out.regenie_log_parser_jar
      )

      regenie_step1_out_ch = REGENIE_STEP1.out.regenie_step1_out
      regenie_step1_parsed_logs_ch = REGENIE_LOG_PARSER_STEP1.out.regenie_step1_parsed_logs

    } else {
      /* 
      You can load pre-made regenie level 1 preds. 
      You must specify a path of /my/path/regenie_step1_out* 
      A file named regenie_step1_out_pred.list must be present together with files like
      regenie_step1_out_1.loco.gz, regenie_step1_out_2.loco.gz, ... (one per phenotype)
      These can be used in STEP 2 only if phenotypes and covariates are exactly the same 
      used to generate step1 predictions (both files and column designation must match exactly)
      */
      regenie_step1_out_ch = Channel
        .fromPath(params.regenie_premade_predictions, checkIfExists: true)
        .collect()

      regenie_step1_parsed_logs_ch = Channel.fromPath("NO_LOG")

    }

    //==== REGENIE STEP 2 ====
    //PARALLELIZE BY CHROM
    if (params.step2_split_by == 'chr') { 
      chromosomes = Channel.of(params.chromosomes)
      bychr_imputed_ch = imputed_plink2_ch.combine(chromosomes)
      REGENIE_STEP2 (
        regenie_step1_out_ch.collect(),
        bychr_imputed_ch,
        VALIDATE_PHENOTYPES.out.phenotypes_file_validated,
        sample_file,
        covariates_file_validated
      )
    //PARALLELIZE BY CHUNK
    } else if (params.step2_split_by == 'chunk') {
      MAKE_CHUNKS(snplist_ch)
      chunks_ch = MAKE_CHUNKS.out.splitText() { it.trim() }
      bychunk_imputed_ch = imputed_plink2_ch.combine(chunks_ch)
      REGENIE_STEP2 (
        regenie_step1_out_ch.collect(),
        bychunk_imputed_ch,
        VALIDATE_PHENOTYPES.out.phenotypes_file_validated,
        sample_file,
        covariates_file_validated
      )
    }

    REGENIE_LOG_PARSER_STEP2 (
      REGENIE_STEP2.out.regenie_step2_out_log.first(),
      CACHE_JBANG_SCRIPTS.out.regenie_log_parser_jar
    )

  //concat by chromosome results into a single result file per pheno
  concat_input_ch = REGENIE_STEP2.out.regenie_step2_out.groupTuple().map{ it -> return tuple(it[0], it[1].flatten())}
  CONCAT_STEP2_RESULTS(concat_input_ch)

  // generate a tuple of phenotype and corresponding result
  CONCAT_STEP2_RESULTS.out.regenie_results_gz
    .flatten()
    .map { it -> return tuple(it.simpleName, it) }
    .set { regenie_step2_by_phenotype }

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
    CLUMP_RESULTS(regenie_step2_by_phenotype, genes_ranges_hg19, genes_ranges_hg38, imputed_plink2_ch, sample_file)
    clump_results_ch = CLUMP_RESULTS.out.best_loci
  } else {
    clump_results_ch = regenie_step2_by_phenotype.map { it -> return tuple(it[0], file('NO_CLUMP_FILE'))}
  }

  merged_results_and_annotated_filtered = regenie_step2_by_phenotype
    .join(ANNOTATE_FILTERED.out.annotated_ch, by: 0)
    .join(clump_results_ch, by: 0, remainder: true)

  //==== GENERATE REPORTS ====
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

workflow.onComplete {
  //==== SAVE CONFIGURATION ====
  pipeline_log_dir = file("$outdir/analysis_config")
  pipeline_log_dir.mkdirs()
  
  def msg="""\
  ## Pipeline information ##
  Pipeline version = ${workflow.manifest.version}
  Pipeline repo = ${workflow.manifest.homePage}
  Nextflow version required = ${workflow.manifest.nextflowVersion}

  ## Execution summary ##
  Execution status: ${ workflow.success ? 'OK' : 'failed' }
  Execution name: ${workflow.runName}
  Execution ID: ${workflow.sessionId}
  Execution date: ${workflow.start}
  Executed by: ${workflow.userName}
  Commandline: ${workflow.commandLine} 
  Container: ${workflow.container}
  """
  .stripIndent()

  params_log = file("${pipeline_log_dir}/pipeline_execution.log")
  params_log.append(msg)

  for (f in workflow.configFiles) {
    config_file = file(f)
    f.copyTo("$pipeline_log_dir")
  }

  nextflow_log = file("${workflow.launchDir}/.nextflow.log")
  nextflow_log.copyTo("${pipeline_log_dir}/nextflow.log")

  //if (params.master_log_dir != null) {
    master_log = file("${master_log_dir}/job_execution_summary.log")
    def master_log_msg="${params.project}\t${workflow.launchDir}\t${ workflow.success ? 'OK' : 'failed' }\n"
    master_log.append(master_log_msg)
  //}

  //CLOSE MESSAGE
  println "Results location: ${ outdir }"
}

workflow.onError {
  if (params.master_log_dir != null) {
    error_log = file("${master_log_dir}/${project_id}_error.log")
    def error_log_msg="""\
    ### ERROR MESSAGE ###
    ${workflow.errorMessage}

    ### ERROR REPORT ###
    ${workflow.errorReport}
    """
    .stripIndent()
    error_log.append(error_log_msg)

    master_log = file("${master_log_dir}/job_execution_summary.log")
    def master_log_msg="${params.project}\t${workflow.launchDir}\t${ workflow.success ? 'OK' : 'failed' }\n"
    master_log.append(master_log_msg)
  }
}