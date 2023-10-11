include { CONVERT_TO_BED        } from '../modules/local/variant_clumping'
include { PLINK_CLUMPING        } from '../modules/local/variant_clumping'
include { MERGE_CLUMP_RESULTS   } from '../modules/local/variant_clumping'
include { MERGE_BED_DATASET     } from '../modules/local/variant_clumping'
include { CHECK_CHANNEL_SIZE    } from '../modules/local/check_channel_size'

workflow CLUMP_RESULTS {
    take:
        results //[val(project_id), val(phenotype), file(pheno_results_gz)]
        genes_interval_hg19 //file
        genes_interval_hg38 //file
        processed_gwas_genotypes //[val(filename), file(bed_bgen_pgen), file(bim_bgi_pvam), file(fam_sample_psam), val(chrom)]

    main:
        //When no LD panel provided we use directly the imputed data to perform clumping
        if (params.ld_panel == 'NO_LD_FILE') {
            //remap to [chrom, file(bed_bgen_pgen), file(bim_bgi_pvam), file(fam_sample_psam)]
            input_ch = processed_gwas_genotypes
                .map { tuple(it[4], it[1], it[2], it[3]) }

            //convert to bed if needed to speed up plink1.9 clumping
            if (params.genotypes_imputed_format != 'bed') {   
                CONVERT_TO_BED(input_ch)
                bed_input_ch = CONVERT_TO_BED.out
            } else {
                bed_input_ch = input_ch
            }   
            
            chromosomes_ch = Channel.fromList(params.chromosomes)
            extension = params.input_format == 'vcf' ? 'vcf.gz' : "${params.input_format}"
            filename = params.input_files.replaceFirst(/(\.bgen|\.vcf\.gz|\.bcf)$/, "")
            
            //Count number of input files
            input_files = []
            input_files = input_files + file("${filename}.${extension}").name
            
            //When input is by chromosome we can use the converted BED directly
            if (params.input_files =~ /\{CHROM\}/) {
                ld_panel_files_ch = bed_input_ch
            } else if (input_files.size() > 1) { 
            //Otherwise, if we have more than 1 input file we need to merge converted BED
                merge_bed_files = bed_input_ch.map{ it[1] }.collect()
                merge_bim_files = bed_input_ch.map{ it[2] }.collect()
                merge_fam_files = bed_input_ch.map{ it[3] }.collect()

                MERGE_BED_DATASET(merge_bed_files, merge_bim_files, merge_fam_files)
                
                merged_bed_by_chrom = chromosomes_ch.combine(MERGE_BED_DATASET.out)
                
                //mix the channels to ensure we have all files
                ld_panel_files_ch = merged_bed_by_chrom 
            } else {
            //Finally, if input is not by chromosome and we have only 1 input file we can use it directly
                ld_panel_files_ch = chromosomes_ch.combine(bed_input_ch.map { tuple(it[1], it[2], it[3]) })
            }          

        } else {
            //If LD panel provided we use it to perform clumping
            //LD panel is expected by chromosome with {CHROM} placeholder
            def pattern = "${params.ld_panel.replace('{CHROM}', '(.+)').replace('/', '\\/')}"
            ld_panel_files_ch = Channel.fromFilePairs("${params.ld_panel.replace('{CHROM}','*')}.{bed,bim,fam}", size:3, flat: true)
                .map { tuple((("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]).replace('.bed',''), it[1], it[2], it[3]) }
        }

        //Filter results channel to keep only chromosomes of interest and check we have all chromosomes
        ld_panel_ch = ld_panel_files_ch.filter { it[0] in params.chromosomes }
        CHECK_CHANNEL_SIZE(ld_panel_ch.count(), params.chromosomes.size(), "LD panel files")        

        //Combine results cahnnel with ld_panel channel
        //Clumping is performed for each chromosome separately to increase speed
        clump_input_ch = results.combine(ld_panel_ch)
        PLINK_CLUMPING(clump_input_ch, genes_interval_hg19, genes_interval_hg38)

        //Now merge all clumpign results by project_id and phenotype
        merge_input_ch = PLINK_CLUMPING.out.chrom_clump_results.groupTuple(by: [0,1], size: params.chromosomes.size())
        MERGE_CLUMP_RESULTS(merge_input_ch)
    
    emit:
        best_loci = MERGE_CLUMP_RESULTS.out.annotloci
}

