import pandas as pd
from glob import glob
import os

files = glob('04_output/wgbs_paired/bismark_meths/cx_report/*.CX_report.txt.gz')

lambdah = {}
chloro  = {}
for f in files:
    print( os.path.basename(f) )
    allc = pd.read_csv(
        f,
        delimiter="\t",
        names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
        )
    # lambda vector
    lambdah[os.path.basename(f)] = allc\
        .loc[allc['chr'] == "Lambda_NEB"]\
        .loc[ (allc['mC'] >0) ]
    lambdah[os.path.basename(f)].insert( 0, 'filename', os.path.basename(f) )
    # chloroplast
    chloro[os.path.basename(f)] = allc\
        .loc[allc['chr'] == "ChrC"]\
        .loc[ (allc['mC'] >0) ]
    chloro[os.path.basename(f)].insert( 0, 'filename', os.path.basename(f) )

outdir = "05_results/01_which_reads_fail/"
pd.concat(lambdah).to_csv(outdir + "lambda_NEB.csv", index = False)
pd.concat(chloro).to_csv( outdir + "ChrC.csv", index = False)