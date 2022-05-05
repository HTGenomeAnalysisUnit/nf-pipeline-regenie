#!/usr/env python
import csv
import sys

input_file = sys.argv[1] # master table from phenotyper script
config_file = sys.argv[2] # template of config file
base_dir = sys.argv[3] # full path to output folder. One sub-folder per run will be created here

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

with open(config_file) as f:
    conf_template = f.read()

for row in read_tsv(input_file):
    run_config_file = open(f"{row['run_group']}/gwas.conf", "w")
    run_template = conf_template
    run_template = update_conf(run_template, 'outdir', base_dir)
    run_template = update_conf(run_template, 'project', row['run_group'])
    run_template = update_conf(run_template, 'covariates_filename', row['cov_file'])
    run_template = update_conf(run_template, 'phenotypes_filename', row['pheno_file'])
    run_template = update_conf(run_template, 'regenie_test', row['genetic_model'])
    run_template = update_conf(run_template, 'phenotypes_binary_trait', row['trait_type'] == "log")
    with open(row['cov_file']) as f:
        first_line = tokenize(f.readline())
        covars = ",".join(first_line[2:])
        run_template = update_conf(run_template, 'covariates_columns', covars)
    
    with open(row['pheno_file']) as f:
        first_line = tokenize(f.readline())
        phenos = ",".join(first_line[2:])
        run_template = update_conf(run_template, 'phenotypes_columns', phenos)

    run_config_file.write(run_template + "\n")
    run_config_file.write(f"manifest.name = 'run-{row['run_group']}-gwas'\n")
    run_config_file.close()