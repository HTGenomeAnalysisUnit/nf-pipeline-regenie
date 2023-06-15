#!/usr/bin/env python

import pandas as pd
import sys

# Read file from command line argument
filename = sys.argv[1]
df = pd.read_csv(filename, sep='\t', header=None)

# Drop columns 20, 21, 22 and group by column 3 (variant)	
df = df.drop(columns=[20, 21, 22])
df[23] = df.groupby(3)[23].transform(lambda x: ','.join(x))

# Convert back to DataFrame and drop duplicates
result = df.drop_duplicates()

# Save result to file
result.to_csv('closegenes_collapsed.tsv', sep='\t', index=False, na_rep='NA')