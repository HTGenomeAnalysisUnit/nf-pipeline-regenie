#!/usr/bin/env python

# Validate input files for Regenie
import os
import sys
import csv
import argparse

REGENIE_MISSING = "NA"

class RegenieValidateInput:
	def __init__(self, input_file, output_file, file_type):
		self.input_file = input_file
		self.output_file = output_file
		self.file_type = file_type

	def validate(self):
		with open(self.input_file, 'r') as input_fh, open(self.output_file, 'w') as output_fh:
			reader = csv.reader(input_fh, delimiter='\t')
			writer = csv.writer(output_fh, delimiter='\t', lineterminator="\n")

			# Get the header row
			header = next(reader)

			# Check if the file is space-separated
			if len(header) == 1:
				input_fh.seek(0)
				reader = csv.reader(input_fh, delimiter=' ')
				header = next(reader)

				# If it's still not valid, print an error message and exit
				if len(header) == 1:
					print(f"ERROR: Input file '{self.input_file}' must be TAB or SPACE separated.", file=sys.stderr)
					return 1

			# Ensure FID and IID are uppercase
			header[0] = header[0].upper()
			header[1] = header[1].upper()

			# Check if the header is valid
			if header[:2] != ['FID', 'IID']:
				print(f"ERROR: header of file '{self.input_file}' must start with 'FID IID'.", file=sys.stderr)
				return 1

			# Write the header to the output file
			writer.writerow(header)

			# Initialize counters for empty and NA values
			count_empty_values = [0] * len(header)
			count_na_values = [0] * len(header)

			# Iterate over the rows in the input file
			for row in reader:
				# Check if the row has the correct number of columns
				if len(row) != len(header):
					print(f"ERROR: Input file '{self.input_file}' parse error in line {reader.line_num}. Detected columns: {len(row)}. Expected columns: {len(header)}.", file=sys.stderr)
					return 1

				# Replace empty values with NA
				if self.file_type == 'phenotype':
					for i in range(2, len(header)):
						if row[i] == '' or row[i] == '.':
							count_empty_values[i] += 1
							row[i] = REGENIE_MISSING
						elif row[i] == 'NA':
							count_na_values[i] += 1

				# Check for empty values in covariate file
				elif self.file_type == 'covariate':
					for i in range(len(row)):
						if row[i] == '':
							print(f"ERROR: Sample {row[0]} includes an empty value in column {i}.", file=sys.stderr)
							return 1

				# Write the row to the output file
				writer.writerow(row)

			# Write the log file
			log_file = os.path.join(os.path.dirname(self.output_file), os.path.splitext(os.path.basename(self.output_file))[0] + '.log')
			with open(log_file, 'w') as log_fh:
				log_writer = csv.writer(log_fh, delimiter='\t', lineterminator="\n")
				log_writer.writerow(['Name', 'Value'])
				log_writer.writerow(['Samples', reader.line_num])
				for i in range(len(count_empty_values)):
					if count_empty_values[i] != 0:
						log_writer.writerow([f"[Phenotype  {header[i]}] NA-replaced empty values", count_empty_values[i]])
					if count_na_values[i] != 0:
						log_writer.writerow([f"[Phenotype  {header[i]}] NA values", count_na_values[i]])

		return 0

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', help='Input file', required=True)
	parser.add_argument('--output', help='Validated output file', required=True)
	parser.add_argument('--type', help='File type', required=True, choices=['phenotype', 'covariate'])
	args = parser.parse_args()

	regenie = RegenieValidateInput(args.input, args.output, args.type)
	exit_code = regenie.validate()
	sys.exit(exit_code)
