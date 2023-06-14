include { CONVERT_TO_BED        } from '../modules/local/variant_clumping'
include { PLINK_CLUMPING        } from '../modules/local/variant_clumping'
include { MERGE_CLUMP_RESULTS   } from '../modules/local/variant_clumping'
include { CHECK_CHANNEL_SIZE    } from '../modules/local/check_channel_size'

workflow CLUMP_RESULTS {
    take:
        results //tuple val(phenotype), file(pheno_results_gz)
        genes_interval_hg19 //file
        genes_interval_hg38 //file
        processed_gwas_genotypes //[val(filename), file(bed_bgen_pgen), file(bim_bgi_pvam), file(fam_sample_psam), val(chrom)]

    main:
        if (params.ld_panel == 'NO_LD_FILE') {
            if (params.genotypes_imputed_format != 'bed') {
                CONVERT_TO_BED(processed_gwas_genotypes)
                bed_files_ch = CONVERT_TO_BED.out
            } else {
                bed_files_ch = processed_gwas_genotypes
                    .map { tuple(it[4], it[1], it[2], it[3]) }
            }

            bed_files_ch.branch { 
                single_file: it[0] == 'ONE_FILE'
                split_by_chr: true
            }
            .set { bed_files_fork_ch }

            chromosomes_ch = Channel.of(params.chromosomes)
            ld_panel_part1_ch = chromosomes_ch.combine(bed_files_fork_ch.single_file)
                .map { tuple(it[0],it[1],it[2],it[3]) }
            ld_panel_ch = ld_panel_part1_ch.mix(bed_files_fork_ch.split_by_chr)
        } else {
            def pattern = "${params.ld_panel.replace('{CHROM}', '(.+)').replace('/', '\\/')}"
            ld_panel_ch = Channel.fromFilePairs("${params.ld_panel.replace('{CHROM}','*')}.{bed,bim,fam}", size:3, flat: true)
                .map { tuple((("${it[1]}" =~ /${pattern}/)[ 0 ][ 1 ]).replace('.bed',''), it[1], it[2], it[3]) }
        }
        
        ld_panel_filtered_ch = ld_panel_ch.filter { it[0] in params.chromosomes }

        CHECK_CHANNEL_SIZE(ld_panel_filtered_ch.count(), params.chromosomes.size(), "LD panel files")
        clump_input_ch = results.combine(ld_panel_filtered_ch)
        PLINK_CLUMPING(clump_input_ch, genes_interval_hg19, genes_interval_hg38)
        merge_input_ch = PLINK_CLUMPING.out.chrom_clump_results.groupTuple()
        MERGE_CLUMP_RESULTS(merge_input_ch)
    
    emit:
        best_loci = MERGE_CLUMP_RESULTS.out.annotloci
}

