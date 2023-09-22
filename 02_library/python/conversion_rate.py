import epiclinestools as epi
import argparse

# Parameters
parser = argparse.ArgumentParser()
parser.add_argument('--input', help="""
Path to the input cytosine report file.
This should be the output of a Bismark cytosine (CX) report, with a row for each cytosine and columns for chromosome, positions, strand, number of methylated reads, total reads, sequence context, and trinucleotide context.
See the help page for `bismark_methylation_extracter` for more, specifically the option `--CX_report`."""
)
parser.add_argument('--output', help="""
File to save the output."""
)
args = parser.parse_args()

cx_report = epi.CytosineCoverageFile(args.input)
conversion = cx_report.conversion_rate()
conversion.to_csv(args.output, index = False)