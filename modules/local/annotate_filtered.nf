process ANNOTATE_FILTERED {

  publishDir "${params.outdir}/results/tophits", mode: 'copy'

  input:
    tuple val(phenotype), path(regenie_merged)
    path genes_hg19, stageAs: 'hg19_genes.bed'
    path genes_hg38, stageAs: 'hg38_genes.bed'

  output:
    tuple val(phenotype), path("${regenie_merged.baseName}.annotated.txt.gz"), emit: annotated_ch

  script:
  def genes = params.genotypes_build == 'hg19' ? "hg19_genes.bed" : "hg38_genes.bed"

  """
  #!/bin/bash
  set -e
  mkdir -p work
  #save original header
  zcat ${regenie_merged} | head -1 > header.txt
  #sort and transform to bed file
  zcat ${regenie_merged} | tail -n+2 | sort -T work -k1,1V -k2,2n | awk '{\$2 = \$2-1 OFS \$2} 1' OFS='\t' > ${regenie_merged.baseName}.sorted.bed
  # bedtools closest will crash if chr names are not the same between -a and -b files
  # we need to provide a genome file
  cut -f1 $genes | uniq | awk '{OFS="\t"}; {print \$1,"1"}' > genome.txt
  # annotate closest gene with bedtools
  bedtools closest -a ${regenie_merged.baseName}.sorted.bed -b ${genes} -d -g genome.txt > ${regenie_merged.baseName}.annotated.bed
  rm ${regenie_merged.baseName}.sorted.bed
  # generate an interval around SNPs to annotate genes
  awk '{OFS="\t"};{\$4=\$3"_%SEP%_"\$4; \$2=\$2-${params.annotation_interval_kb * 1000}; \$3=\$3+${params.annotation_interval_kb * 1000}}; {print \$0}' ${regenie_merged.baseName}.annotated.bed \
    | awk '{OFS="\t"}; \$2 < 0 {\$2 = 0}; {print ;}' \
    | bedtools intersect -a stdin -b $genes -loj \
    | cut -f1,4- | sed 's/_%SEP%_/\t/' \
    | awk '{\$2 = \$2-1 OFS \$2} 1' OFS='\t' \
    > ${regenie_merged.baseName}.closegenes.bed
  # merge resulting bed to get annotation by SNP
  bedtools merge -i ${regenie_merged.baseName}.closegenes.bed -c \$(echo {4..20} | tr " " ","),24 -o \$(printf 'distinct,%.0s' {1..18}) > ${regenie_merged.baseName}.final.bed
  rm ${regenie_merged.baseName}.closegenes.bed
  # remove duplication of 2nd column
  cut -f1,3- ${regenie_merged.baseName}.final.bed > ${regenie_merged.baseName}.final.fixed.bed
  rm ${regenie_merged.baseName}.final.bed
  # write extended header
  (cat header.txt | sed ' 1 s/.*/&\tCLOSEST_GENE_CHROMOSOME\tCLOSEST_GENE_START\tCLOSEST_GENE_END\tCLOSEST_GENE_NAME\tCLOSEST_GENE_DISTANCE\tGENES_${params.annotation_interval_kb}KB/' && cat ${regenie_merged.baseName}.final.fixed.bed) > ${regenie_merged.baseName}.final.header.bed
  rm ${regenie_merged.baseName}.final.fixed.bed
  # sort by p-value again
  (cat ${regenie_merged.baseName}.final.header.bed | head -n 1 && cat ${regenie_merged.baseName}.final.header.bed | tail -n +2 | sort -T work -k13,13gr) | gzip > ${regenie_merged.baseName}.annotated.txt.gz
  """
}
