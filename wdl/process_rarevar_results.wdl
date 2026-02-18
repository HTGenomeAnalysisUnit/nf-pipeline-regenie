version 1.0

## ============================================================================
## Subworkflow: Process Rare Variant Results (filter + corrected p-values)
## ============================================================================

import "structs.wdl" as Structs
import "tasks.wdl" as Tasks

workflow ProcessRarevarResults {
    input {
        String project_id
        Array[String] phenotypes
        Array[File] merged_results_by_pheno  # parallel with phenotypes
        Float rarevar_min_log10p = 5.0
    }

    scatter (pheno_idx in range(length(phenotypes))) {
        String phenotype = phenotypes[pheno_idx]
        File pheno_results = merged_results_by_pheno[pheno_idx]

        # Filter results
        call Tasks.FilterResults as FilterRarevarResults {
            input:
                project_id = project_id,
                phenotype = phenotype,
                regenie_result_gz = pheno_results,
                annotation_min_log10p = rarevar_min_log10p,
                rarevar_results = true
        }

        # Process rare variant results (corrected p-values)
        call Tasks.ProcessRarevarResults {
            input:
                project_id = project_id,
                phenotype = phenotype,
                regenie_result_gz = pheno_results
        }
    }

    output {
        Array[File] filtered_results = FilterRarevarResults.filtered_results
        Array[File] processed_results = ProcessRarevarResults.processed_results
    }
}
