#!/usr/bin/env python
import csv
import sys
import os

input_file = sys.argv[1] # master table from phenotyper script
#config_file = sys.argv[2] # template of config file
#base_dir = sys.argv[3] # full path to output folder. One sub-folder per run will be created here

def read_tsv(tsv_file):
    """
    Reads a tsv file and returns a list of dictionaries.
    """
    with open(tsv_file, 'r') as tsv_in:
        reader = csv.DictReader(tsv_in, delimiter='\t')
        return [row for row in reader]

def tokenize(x, sep="\t"):
    out = x.strip("\n")
    out = out.split(sep)
    return out

def update_conf(template, tag, value):
    if isinstance(value, bool):
        config_value = 'true' if value else 'false'
    else:
        config_value = f"'{value}'"

    template = template.replace(tag, f"{tag} = {config_value}")
    return template

def find_not_overlap(lst1, lst2):
    lst3 = [value for value in lst1 if value not in lst2]
    return lst3

# with open(config_file) as f:
#     conf_template = f.read()

#subfolder with tables has the same name as run_group
#then in this folder we have chunk_<i>.cov chunk_<i>.pheno
#Read master table and output a new table with columns reflecting the data channel we need

with open("analysis.conf", "w") as outf:
    outf.write('\t'.join(['project_id','pheno_file','pheno_cols','pheno_binary','pheno_model','cov_file','cov_cols','cov_cat_cols']) + '\n')
    for row in read_tsv(input_file):
        if row['cat_var'] != 'NA':
            cat_covars = row['cat_var'].split(",")
        else:
            cat_covars = []

        covars = ''
        covars_file = 'NO_COV_FILE'
        if row['cov_file'] != "NO_COV_FILE":
            covars_file = os.path.abspath(row['cov_file'])
            with open(row['cov_file']) as f:
                first_line = tokenize(f.readline())
                covars = first_line[2:]
                covars = find_not_overlap(covars, cat_covars)
                covars = ",".join(covars)

        with open(row['pheno_file']) as f:
            first_line = tokenize(f.readline())
            phenos = ",".join(first_line[2:])
                
        outf.write('\t'.join([
            row['run_group'],
            os.path.abspath(row['pheno_file']),
            phenos,
            str(row['trait_type'] == "log"),
            row['genetic_model'],
            covars_file,
            covars,
            row['cat_var']]) + '\n')