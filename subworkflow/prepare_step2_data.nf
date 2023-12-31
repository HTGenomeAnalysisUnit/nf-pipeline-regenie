include { CONVERT_TO_BGEN   } from '../modules/local/imputed_to_plink2' addParams(outdir: "${params.outdir}/converted_bgen", publish: params.save_converted_bgen, dosage_from: params.dosage_from)
include { MAKE_BGEN_INDEX   } from '../modules/local/make_bgen_index'   addParams(outdir: "${params.outdir}/bgen_dataset", publish: params.save_bgen_index)
include { MAKE_BGEN_SAMPLE  } from '../modules/local/make_bgen_sample'  addParams(outdir: "${params.outdir}/bgen_dataset", publish: params.save_bgen_sample)
include { CHECK_MAX_CHANNEL_SIZE } from '../modules/local/check_channel_size'

workflow PREPARE_GENETIC_DATA {
    main:
    //First we check in which way input data is provided 
    //In the end we want to have a channel of tuples 
    //[val(basename), file(bed/bgen/pgen), file(bim/bgi/pvar), file(fam/sample/psam), val(chromosome/one_file)]
    
    //==== INPUT DATA IS PROVIDED SPLITTED BY CHROMOSOMES ====
    if (params.genotypes_data =~ /\{CHROM\}/) {
        def pattern = "${params.genotypes_data.replace('{CHROM}', '(.+)').replace('/', '\\/')}"
        switch(params.input_format) {
        case "bcf":
            input_ch = Channel.fromFilePairs("${params.genotypes_data.replace('{CHROM}','*')}", size: 1, flat: true)
                .map { tuple(it[0], it[1], (("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]).replace('.bcf','')) }
            genotypes_files = input_ch.filter { it[2] in params.chromosomes }
            break
        case "vcf":
            input_ch = Channel.fromFilePairs("${params.genotypes_data.replace('{CHROM}','*')}", size: 1, flat: true)
                .map { tuple(it[0], it[1], (("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]).replace('.vcf.gz','')) }
            genotypes_files = input_ch.filter { it[2] in params.chromosomes }
            break
        case "bgen":
            input_ch = Channel.fromFilePairs("${params.genotypes_data.replace('{CHROM}','*')}", size: 1, flat: true)
                .map { tuple(it[0], it[1], file("${it[1]}.bgi"), file("${it[1].parent}/${it[1].baseName}.sample"), (("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]).replace('.bgen','')) }
            genotypes_files = input_ch.filter { it[4] in params.chromosomes }
            break
        case "bed":
            input_ch = Channel.fromFilePairs("${params.genotypes_data.replace('{CHROM}','*')}.{bed,bim,fam}", size:3, flat: true)
                .map { tuple(it[0], it[1], it[2], it[3], (("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]).replace('.bed','')) }
            genotypes_files = input_ch.filter { it[4] in params.chromosomes }
            break
        case "pgen":
            input_ch = Channel.fromFilePairs("${params.genotypes_data.replace('{CHROM}','*')}.{pgen,psam,pvar}", size:3, flat: true)
                .map { tuple(it[0], it[1], it[3], it[2], (("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]).replace('.pgen','')) }
            genotypes_files = input_ch.filter { it[4] in params.chromosomes }
            break
        default:
            log.error "Unknown input format: ${params.input_format}"
            exit 1
        } 
    
    } else {
    //==== INPUT DATA IS PROVIDED IN A SINGLE FILE ==== 
        switch(params.input_format) {
        case "bcf":
            genotypes_files = Channel.fromFilePairs("${params.genotypes_data}", checkIfExists: true, size: 1, flat: true)
                .map { tuple(it[0], it[1], "ONE_FILE") }
            break
        case "vcf":
            genotypes_files = Channel.fromFilePairs("${params.genotypes_data}", checkIfExists: true, size: 1, flat: true)
                .map { tuple(it[0], it[1], "ONE_FILE") }
            break
        case "bgen":
            genotypes_files = Channel.fromFilePairs("${params.genotypes_data}", checkIfExists: true, size: 1, flat: true)
                .map { tuple(it[0], it[1], file("${it[1]}.bgi"), file("${it[1].parent}/${it[1].baseName}.sample"), "ONE_FILE") }
            break
        case "bed":
            genotypes_files = Channel.fromFilePairs("${params.genotypes_data}.{bed,bim,fam}", checkIfExists: true, size:3, flat: true)
                .map { tuple(it[0], it[1], it[2], it[3], "ONE_FILE") }
            break
        case "pgen":
            genotypes_files = Channel.fromFilePairs("${params.genotypes_data}.{pgen,psam,pvar}", checkIfExists: true, size:3, flat: true)
                .map { tuple(it[0], it[1], it[3], it[2], "ONE_FILE") }
            break
        default:
            log.error "Unknown input format: ${params.input_format}"
            exit 1
        }
        //Check we have a single input file
        //CHECK_MAX_CHANNEL_SIZE(genotypes_files.count(), 1, "genotypes_files")
    }
    
    //==== SET OUTPUT CHANNEL ====
    genotypes_plink2_ch = genotypes_files

    //==== CONVERT TO PGEN IF INPUT IS VCF ====
    if (params.input_format in ["vcf","bcf"]) {
        CONVERT_TO_BGEN ( genotypes_files )
        genotypes_plink2_ch = CONVERT_TO_BGEN.out.genotypes_data
    } 
    
    if (params.input_format == "bgen") {
        //==== MAKE BGI INDEX IF MISSING ====
        genotypes_files.map { tuple(it[0], it[1], it[2], it[4]) }
            .branch {
                    found: it[2].exists()
                    missing: true
                }
            .set { check_bgi_ch }
        
        MAKE_BGEN_INDEX(check_bgi_ch.missing.map { tuple(it[0], it[1], it[3]) } )
        bgen_bgi_ch = check_bgi_ch.found
            .mix(MAKE_BGEN_INDEX.out)

        //==== CHECK SAMPLE FILE FOR BGEN ====
        //If a sample file is provided in params, this is used for all bgen files
        if (params.bgen_sample_file != 'NO_SAMPLE_FILE') {
            genotypes_bgen_and_sample = bgen_bgi_ch
                .map{ tuple(it[0], it[1], it[2], file(params.bgen_sample_file, checkIfExists: true), it[3]) }
        } else {
        //Otherwise check that a sample file exists and make one when missing
            genotypes_files.map { tuple(it[0], it[1], it[3], it[4]) }
            .branch {
                    found: it[2].exists()
                    missing: true
                }
            .set { check_sample_ch }

            MAKE_BGEN_SAMPLE(check_sample_ch.missing.map { tuple(it[0], it[1], it[3]) } )
            
            bgen_sample_ch = check_sample_ch.found
                .mix(MAKE_BGEN_SAMPLE.out)

            //put together bgi index and sample file for all datasets
            genotypes_bgen_and_sample = bgen_bgi_ch.join(bgen_sample_ch, by: [0,3])
                .map { tuple(it[0], it[2], it[3], it[5], it[1]) }
        }

        genotypes_plink2_ch = genotypes_bgen_and_sample
    }

    emit:
        processed_genotypes = genotypes_plink2_ch //[file_prefix, file(bed_bgen_pgen), file(bim_bgi_pvar), file(fam_sample_psam), val(chrom)]
}