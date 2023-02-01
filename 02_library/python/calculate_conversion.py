import pandas as pd
import numpy as np
import argparse

# Parameters
parser = argparse.ArgumentParser()
parser.add_argument('--allc', help="""
Path to the input allc file.
This should be the output of a Bismark cytosine report, with a row for each cytosine and columns for chromosome, positions, strand, number of methylated reads, total reads, sequence context, and trinucleotide context.
See the help page for `bismark_methylation_extracter` for more, specifically the option `--cytosine_report`."""
)
parser.add_argument('--output', help="""
File to save the output."""
)
args = parser.parse_args()

# Import the allc file
allc_df = pd.read_csv(
    args.allc,
    delimiter="\t",
    names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
    )

# Conversion rate on each chromosome
def conversion_rate(allc_file):
    return allc_file.groupby('chr', group_keys=False).apply(
        lambda x : float(x['mC'].sum()) / (x['uC'] + x['mC']).sum() 
        )
conversion = conversion_rate(allc_df)

print(conversion)

# Format header and data lines to output
header_names = '\t'.join(conversion.keys())
outline = '\t'.join(
    [str(round(item, 3)) for item in conversion.to_list() ]
    )

# Write to disk
with open(args.output, 'a') as outfile:
    if outfile.tell() == 0:
        outfile.write(header_names)
    outfile.write('\n' + outline)