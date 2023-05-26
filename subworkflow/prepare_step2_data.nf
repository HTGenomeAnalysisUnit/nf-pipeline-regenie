def check_size(file_glob, max_size) {
    def files = file(file_glob)
    def size = files.size()
    if (size > max_size) {
        println "ERROR: ${file_glob} returns ${size} files, but only ${max_size} are allowed"
        exit 1
    }
}

include { CONVERT_TO_PGEN   } from '../modules/local/imputed_to_plink2' addParams(outdir: "${params.outdir}/converted_PGEN", publish: params.save_pgen)
include { MAKE_BGEN_INDEX  } from '../modules/local/make_bgen_index' addParams(outdir: "${params.outdir}/bgi_index", publish: params.save_bgen_index)
include { MAKE_SNPLIST      } from '../modules/local/make_snplist' addParams(outdir: "${params.outdir}/snplist", publish: params.save_snplist)

workflow PREPARE_GENETIC_DATA {
    take:
    chromosomes //a list of chromosome names to process
    
    main:
    //First we check in which way input data is provided 
    //In the end we want to have a channel of tuples 
    //[val(basename), file(bed/bgen/pgen), file(bim/bgi/pvar), file(fam/sample/psam), val(chromosome/one_file)]
    
    //==== INPUT DATA IS PROVIDED SPLITTED BY CHROMOSOMES ====
    if (params.genotypes_data =~ /\{CHROM\}/) {
        def pattern = "${params.genotypes_data.replace('{CHROM}', '(.+)').replace('/', '\\/')}"
        switch(params.input_format) {
        case "vcf":
            input_ch = Channel.fromFilePairs("${params.genotypes_data.replace('{CHROM}','*')}", size: 1, flat: true)
                .map { tuple(it[0], it[1], ("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]) }
            break
        case "bgen":
            input_ch = Channel.fromFilePairs("${params.genotypes_data.replace('{CHROM}','*')}", size: 1, flat: true)
                .map { tuple(it[0], it[1], file("${it[1]}.bgi"), ("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]) }
            break
        case "bed":
            input_ch = Channel.fromFilePairs("${params.genotypes_data.replace('{CHROM}','*')}.{bed,bim,fam}", size:3, flat: true)
                .map { tuple(it[0], it[1], it[2], it[3], ("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]) }
            break
        case "pgen":
            input_ch = Channel.fromFilePairs("${params.genotypes_data.replace('{CHROM}','*')}.{pgen,pvar,psam}", size:3, flat: true)
                .map { tuple(it[0], it[1], it[2], it[3], ("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]) }
            break
        default:
            println "Unknown input format: ${params.input_format}"
            exit 1
        } 

        //Now we filter out chromosomes that are not in the list
        genotypes_files = input_ch.filter { it[2] in chromosomes }
    
    } else {
    //==== INPUT DATA IS PROVIDED IN A SINGLE FILE ====
        //Check we have a single input file
        check_size(params.genotypes_data, 1)
       
        switch(params.input_format) {
        case "vcf":
            genotypes_files = Channel.fromFilePairs("${params.genotypes_data}", checkIfExists: true, size: 1, flat: true)
                .map { tuple(it[0], it[1], "ONE_FILE") }
            break
        case "bgen":
            genotypes_files = Channel.fromFilePairs("${params.genotypes_data}", checkIfExists: true, size: 1, flat: true)
                .map { tuple(it[0], it[1], file("${it[1]}.bgi"), "ONE_FILE") }
            break
        case "bed":
            genotypes_files = Channel.fromFilePairs("${params.genotypes_data}.{bed,bim,fam}", checkIfExists: true, size:3, flat: true)
                .map { tuple(it[0], it[1], it[2], it[3], "ONE_FILE") }
            break
        case "pgen":
            genotypes_files = Channel.fromFilePairs("${params.genotypes_data}.{pgen,pvar,psam}", checkIfExists: true, size:3, flat: true)
                .map { tuple(it[0], it[1], it[2], it[3], "ONE_FILE") }
            break
        default:
            println "Unknown input format: ${params.input_format}"
            exit 1
        }
    }
    
    //==== CONVERT TO PGEN IF INPUT IS VCF ====
    if (params.input_format == "vcf") {
        CONVERT_TO_PGEN ( genotypes_files )
        genotypes_plink2_ch = CONVERT_TO_PGEN.out.genotypes_data
    } else {
        genotypes_plink2_ch = genotypes_files
    }

    
    if (params.input_format == "bgen") {
        //==== MAKE BGI INDEX IF MISSING ====
        genotypes_files
            .branch {
                    found: it[2].exists()
                    missing: true
                }
            .set { check_bgi_ch }
        
        MAKE_BGEN_INDEX(check_bgi_ch.missing.map { tuple(it[0], it[1], it[3]) } )
        bgen_bgi_ch = check_bgi_ch.found
            .mix(MAKE_BGEN_INDEX.out.bgen_with_index)

        //==== CHECK SAMPLE FILE FOR BGEN ====
        file("${it[1].parent}/${it[1].baseName}.sample")
        //If a sample file is provided in params, this is used for all bgen files
        if (params.bgen_sample_file != 'NO_SAMPLE_FILE') {
            genotypes_bgen_and_sample = bgen_bgi_ch
                .map{ tuple(it[0], it[1], it[2], file(params.bgen_sample_file, checkIfExists: true), it[3]) }
        } else {
        //Otherwise check that a sample file exists and make one when missing
            bgen_bgi_ch
            .map { tuple(it[0], it[1], it[2], file("${it[1].parent}/${it[1].baseName}.sample"), it[3]) }
            .branch {
                    found: it[3].exists()
                    missing: true
                }
            .set { check_sample_ch }   
        }

        MAKE_BGEN_SAMPLE(check_sample_ch.missing.map { tuple(it[0], it[1], it[2], it[4]) } )
        genotypes_plink2_ch = check_sample_ch.found
            .mix(MAKE_BGEN_SAMPLE.out.bgen_with_sample)
    }

    emit:
        processed_genotypes = genotypes_plink2_ch
}