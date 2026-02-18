version 1.0

## ============================================================================
## Subworkflow: Process GWAS Results (filter, annotate, optional clumping)
## ============================================================================

import "structs.wdl" as Structs
import "tasks.wdl" as Tasks

workflow ProcessGwasResults {
    input {
        String project_id
        Array[String] phenotypes
        Array[File] merged_results_by_pheno  # parallel with phenotypes

        # Gene annotation files
        File genes_bed
        File? genes_ranges  # for clumping

        # Annotation params
        Float annotation_min_log10p = 7.3
        Int annotation_interval_kb = 25
        File collapse_closegenes_py

        # Clumping params
        Boolean clumping = true
        Float clump_p1 = 5e-8
        Float clump_p2 = 1e-4
        Int clump_kb = 250
        String genotypes_build
        String genotypes_imputed_format = "bgen"

        # LD panel files (one per chromosome, or use imputed data)
        Array[String] chromosomes_list
        String ld_panel = "NO_LD_FILE"
        Array[File]? ld_panel_beds
        Array[File]? ld_panel_bims
        Array[File]? ld_panel_fams
        Array[String]? ld_panel_chroms

        # Imputed genotype data (for clumping when no LD panel)
        Array[GenotypeFiles]? processed_gwas_genotypes
    }

    # Process each phenotype
    scatter (pheno_idx in range(length(phenotypes))) {
        String phenotype = phenotypes[pheno_idx]
        File pheno_results = merged_results_by_pheno[pheno_idx]

        # Filter results
        call Tasks.FilterResults as FilterGwasResults {
            input:
                project_id = project_id,
                phenotype = phenotype,
                regenie_result_gz = pheno_results,
                annotation_min_log10p = annotation_min_log10p,
                rarevar_results = false
        }

        # Annotate filtered results
        call Tasks.AnnotateFiltered {
            input:
                project_id = project_id,
                phenotype = phenotype,
                regenie_merged = FilterGwasResults.filtered_results,
                genes_bed = genes_bed,
                annotation_interval_kb = annotation_interval_kb,
                collapse_closegenes_py = collapse_closegenes_py
        }

        # Clumping per phenotype (if enabled)
        if (clumping && defined(genes_ranges)) {
            # When LD panel is provided, scatter over chromosomes
            if (ld_panel != "NO_LD_FILE" && defined(ld_panel_beds)) {
                scatter (chrom_idx in range(length(select_first([ld_panel_chroms, []])))) {
                    call Tasks.PlinkClumping as ClumpWithLdPanel {
                        input:
                            project_id = project_id,
                            phenotype = phenotype,
                            pheno_results_gz = pheno_results,
                            chromosome = select_first([ld_panel_chroms, []])[chrom_idx],
                            bed = select_first([ld_panel_beds, []])[chrom_idx],
                            bim = select_first([ld_panel_bims, []])[chrom_idx],
                            fam = select_first([ld_panel_fams, []])[chrom_idx],
                            genes_interval = select_first([genes_ranges]),
                            genotypes_build = genotypes_build,
                            clump_p1 = clump_p1,
                            clump_p2 = clump_p2,
                            clump_kb = clump_kb,
                            annotation_interval_kb = annotation_interval_kb
                    }
                }

                call Tasks.MergeClumpResults as MergeClumps {
                    input:
                        project_id = project_id,
                        phenotype = phenotype,
                        chromosome_clumps = select_all(ClumpWithLdPanel.clumped),
                        chromosome_ranges = select_all(ClumpWithLdPanel.clumped_ranges)
                }
            }
        }
    }

    output {
        Array[File] filtered_results = FilterGwasResults.filtered_results
        Array[File] annotated_results = AnnotateFiltered.annotated_results
        Array[File?] toploci = MergeClumps.toploci
        Array[File?] annotloci = MergeClumps.annotloci
    }
}
