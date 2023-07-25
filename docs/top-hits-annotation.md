## Annotation of tophits

For top hits SNPs / loci the pipeline performs automatic annotation with overlapping genes and nearby genes. 

The maximum distance between the variant / locus and the sorrounding annotated gene is set with the parameter `annotation_interval_kb` (25kb by default). Genes located +/- defined kb from the SNP / locus are annotated as nearby genes.

Pre-made BED files for gene annotation are available for GRCh37 and GRCh38. These files are generated from the [GENCODE annotation](https://www.gencodegenes.org/) version 39. When performing annotation one can chose to annotate for all genes or only protein coding genes by setting the `genes_group` parameter to `all` or `protein_coding` respectively. By default, the pipeline will annotate for protein coding genes.

If you want to use a different annotation file, you can provide custom gene definition with the `genes_bed` and `genes_ranges` parameters. In this case the pre-made annotation files will not be used. For this you need to provide:

- `genes_bed`: a regular tab-separated BED file (chromosome, start, end; zero-based), with gene names in the 4th column. Note that a name must be present for all genes.
- `genes_ranges`: essentially the same as above but space separated.