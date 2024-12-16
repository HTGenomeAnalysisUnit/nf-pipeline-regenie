#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
from statsmodels.sandbox.stats.multicomp import multipletests

sumstat_file = sys.argv[1]

#Compute bonferroni and FDR corrected P across all tests
dt = pd.read_csv(sumstat_file, sep='\t', skiprows=1)
dt['ALLELE1'] = dt['ALLELE1'].astype(str)
dt['LOG10P_BONF_alltests'] = -np.log10(multipletests(10 ** -dt['LOG10P'], method='bonferroni')[1])
dt['LOG10P_FDR_alltests'] = -np.log10(multipletests(10 ** -dt['LOG10P'], method='fdr_bh')[1])

# Compute bonferroni and FDR corrected P for each specific test
# This means each combination of AF bin, mask and statistical test
dt['LOG10P_BONF_bygroup'] = dt.groupby(['ALLELE1', 'TEST'])['LOG10P'].transform(lambda x: -np.log10(multipletests(10 ** -x, method='bonferroni')[1]))
dt['LOG10P_FDR_bygroup'] = dt.groupby(['ALLELE1', 'TEST'])['LOG10P'].transform(lambda x: -np.log10(multipletests(10 ** -x, method='fdr_bh')[1]))

# Create 2 tables, one for aggregated burden tests and one with others. 
# Any test with TEST containing ADD-BURDEN is a burden test
burden_tests = dt[dt['TEST'].str.contains('ADD-BURDEN')]
other_tests = dt[~dt['TEST'].str.contains('ADD-BURDEN')]

# Split values in the ALLELE1 column into two columns named MASK and AFBIN using the first dot as separator
dt_clean = other_tests.assign(tmp_col=dt['ALLELE1'].apply(lambda x: x.replace('.', '@', 1))) \
  .assign(MASK=lambda x: x['tmp_col'].apply(lambda x: x.split('@')[0])) \
  .assign(AFBIN=lambda x: x['tmp_col'].apply(lambda x: x.split('@')[1])) \
  .drop(columns=['tmp_col'])

# For aggregated tests set MASK and AFBIN to AGGREGATED-BURDEN
burden_tests_clean = burden_tests.assign(MASK='AGGREGATED-BURDEN', AFBIN='AGGREGATED-BURDEN')

# Concatenate the two tables
dt_clean = pd.concat([dt_clean, burden_tests_clean])

# Create a GENE column by taking the first part of the ID column
dt_clean['GENE'] = dt_clean.apply(lambda x: x['ID'].split('.')[0], axis=1)

# Save table with corrected P values
prefix = sumstat_file.rstrip('.gz')
dt_clean.to_csv(f"{prefix}.correctedP.gz", sep='\t', index=False, na_rep='NA')