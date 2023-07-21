#!/usr/bin/env python

# Parse regenie log files into a table
import argparse
import csv
import re

REGENIE_CALL_PATTERN = "Options in effect:"

def main():
	parser = argparse.ArgumentParser(description="Regenie log file")
	parser.add_argument("files", metavar="files", type=str, nargs="+", help="Regenie log file")
	parser.add_argument("--output", metavar="output", type=str, required=True, help="Output file")
	args = parser.parse_args()

	assert args.files is not None
	assert args.output is not None

	with open(args.output, 'w') as output_fh:
		writer = csv.writer(output_fh, delimiter="\t", lineterminator="\n")

		columnsWrite = ["Name", "Value"]
		writer.writerow(columnsWrite)

		warningMsgs = []
		optionsInEffect = []

		countVariants = 0
		countFiles = 0

		for file in args.files:
			countFiles += 1
			isRegenieCall = False

			with open(file, "r") as f:
				for line in f:
					line = line.strip()

					if REGENIE_CALL_PATTERN in line:
						isRegenieCall = True
						continue

					if "-summary :" in line:
						pattern = r'([0-9]+) variants'
						match = re.search(pattern, line)
						if match:
							value = match.group(1).strip()
							countVariants += int(value)
					elif "n_snps =" in line and "pvar" in line:
						value = line.split("=")[1].strip()
						countVariants += int(value)
					elif "WARNING:" in line:
						warningMsgs.append(line)

					if countFiles > 1:
						continue

					if isRegenieCall:
						if not line.endswith("\\"):
							isRegenieCall = False
						else:
							optionsInEffect.append(line[:line.index("\\")])

					if "REGENIE v" in line:
						value = re.split(r'\s+', line)[2].strip()
						writer.writerow(["Regenie Version", value])
					elif "n_snps = " in line and "bim" in line:
						value = line.split("=")[1].strip()
						writer.writerow(["Variants total (*.bim)", int(value)])
					elif "+number of variants remaining in the analysis" in line:
						value = line.split("=")[1].strip()
						writer.writerow(["Variants used (*.bim)", int(value)])
					elif "+number of genotyped individuals to keep" in line:
						value = line.split("=")[1].strip()
						writer.writerow(["Samples total (*.fam)", int(value)])
					elif "n_pheno = " in line:
						value = line.split("=")[1].strip()
						writer.writerow(["Number of defined phenotypes", int(value)])
					elif "n_cov = " in line:
						value = line.split("=")[1].strip()
						writer.writerow(["Number of defined covariates", int(value)])
					elif "-number of phenotyped individuals" in line:
						value = line.split("=")[1].strip()
						writer.writerow(["Phenotyped individuals total", int(value)])
					elif "number of individuals used in analysis" in line:
						value = line.split("=")[1].strip()
						writer.writerow(["Phenotyped individuals used", int(value)])
					elif "--minMAC" in line:
						pattern = r'--minMAC\s+([0-9]+)'
						match = re.search(pattern, line)
						if match:
							value = match.group(1).strip()
							writer.writerow(["MAC limit", value])
					elif "--minINFO" in line and "is skipped" not in line:
						pattern = r'--minINFO\s+([0-9.]+)'
						match = re.search(pattern, line)
						if match:
							value = match.group(1).strip()
							writer.writerow(["Imputation info score limit", value])
					elif "Number of ignored SNPs due to low MAC or info score" in line:
						value = line.split(":")[1].strip()
						writer.writerow(["Variants ignored (low MAC or low info score)", value])
					elif "Number of ignored tests due to low MAC or info score" in line:
						value = line.split(":")[1].strip()
						writer.writerow(["Variants ignored (low MAC or low info score)", value])

		if countVariants > 0:
			writer.writerow(["Variants used (*.bgen or *.pvar)", countVariants])

		if len(args.files) > 1:
			writer.writerow(["Number of parsed log files", len(args.files)])

		if warningMsgs:
			for warn in warningMsgs:
				writer.writerow(["Warnings", warn])

		writer.writerow(["Regenie Call", " ".join(optionsInEffect)])

if __name__ == "__main__":
	main()
