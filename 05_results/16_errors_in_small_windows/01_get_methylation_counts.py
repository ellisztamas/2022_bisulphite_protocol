"""
Script to import the TAIR10 annotation file and pull out information methylated
and unmethylated reads in each context for:
1. genes
2. RdDM- and CMT-regulated TEs
3. windows of 1000 base pairs across each chromosome
for the genotyping file on the 30x Col-0 data we ran.
"""

import epiclinestools as epi
import pandas as pd
import argparse

# Parameters
parser = argparse.ArgumentParser(description = 'Parse parameters ')
parser.add_argument('-i', '--input', help = 'A cytosine report file generated by Bismark', required = True)
parser.add_argument('-w', '--windows', help = 'CSV file', required = True)
# parser.add_argument('-f', '--features', help = 'CSV file', required = True)
args = parser.parse_args()

print("Using epiclinestools verion {}".format(epi.__version__))

# Cytosine coverage file
print("Importing coverage file.")
col0 = epi.CytosineCoverageFile(args.input)

# WINDOWS
# Get counts in windows of 150 bp
print("Counting methylated reads in windows")
col0.methylation_in_windows(window_size = 150).to_csv(
    args.windows, index=False
)