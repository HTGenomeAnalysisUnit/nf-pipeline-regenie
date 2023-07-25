# Variant annotation files

When running a rare variant analysis additional variant annotation files are needed and must be provided using the following parameters:

- `rarevar_set_list_file`: a set list file linking variants to genes. Essentially a tab- or space-separated file with 4 columns: the set/gene name, the chromosome and physical position for the set/gene, a comma-separated list of variants included in the set/gene. This file is used to define the variant sets to be tested. The chromosome names must be numbers 1-22 for autosomes.

- `rarevar_anno_file`: variant annotation file. A tab- or space-separated file with 3 columns: variant name, the set/gene name, a single annotation category (for example missense, LoF, ...). Variants not in this file will be assigned to a default "NULL" category. A maximum of 63 annotation categories (+NULL category) is allowed. 

- `rarevar_mask_file`: mask definition file. A tab- or space-separated file with 2 columns: a mask name followed by a comma-separated list of variant categories included in the mask.

More information on the exact foramt of these files is available in the [regenie documentation](https://rgcgithub.github.io/regenie/options/#annotation-input-files).

For each gene defined in the set list file, regenie will apply the configured tests separately for each mask defined in the mask file (thus including only the variant categories described by each mask).