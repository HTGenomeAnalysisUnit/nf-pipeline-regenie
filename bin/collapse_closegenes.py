#!/usr/bin/env python

import pandas as pd
import shutil
import sys

# Read file from command line argument
filename = sys.argv[1]
try:
	df = pd.read_csv(filename, sep='\t', header=None)

	# Drop columns 20, 21, 22 and group by column 3 (variant)	
	last_column = df.columns[-1]
	cols_to_drop = list(range(last_column-3, last_column))
	df = df.drop(columns=cols_to_drop)
	df[last_column] = df.groupby(3)[last_column].transform(lambda x: ','.join(x))

	# Convert back to DataFrame and drop duplicates
	result = df.drop_duplicates()

	# Save result to file
	result.to_csv('closegenes_collapsed.tsv', sep='\t', index=False, na_rep='NA', header=False)

except pd.errors.EmptyDataError:
	print(f"Failed reading from {filename}. The file is empty")
	shutil.copyfile(filename, "closegenes_collapsed.tsv")