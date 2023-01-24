import pandas as pd
import numpy as np


def conversion_rate(allc_file):
    # 1.0*allc_file['mC'].sum() / (allc_file['mC'] + allc_file['uC']).sum()
    return allc_file.groupby('chr').apply(
        lambda x : float(x['mC'].sum()) / (x['uC'] + x['mC']).sum() 
        )

# Trimming 5-prime ends
allc_path1 = '04_output/wgbs_paired/bismark_meths/cx_report/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz'
allc1 = pd.read_csv(
    allc_path1,
    delimiter="\t",
    names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
    )
conversion1 = conversion_rate(allc1)

allc1.loc[allc1['chr'] == "Lambda_NEB"]\
    .loc[ (allc1['mC'] >0) ]

allc = allc1




# Trimming 5 and 3 prime ends
allc_path2 = '04_output/three_prime_15/bismark_meths/cx_report/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz'
allc2 = pd.read_csv(
    allc_path2,
    delimiter="\t",
    names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
    )
conversion2 = conversion_rate(allc2)


# Trimming 5 prime for both read, and 3 prime  for R1
allc_path3 = '04_output/trim_r1/bismark_meths/cx_report/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz'
allc3 = pd.read_csv(
    allc_path3,
    delimiter="\t",
    names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
    )
conversion3 = conversion_rate(allc3)

# Trim 5prime, keep reads >100bp
allc_path4 = '04_output/wgbs_min_100bp/bismark_meths/cx_report/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz'
allc4 = pd.read_csv(
    allc_path4,
    delimiter="\t",
    names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
    )
conversion4 = conversion_rate(allc4)

# Trim 5 ad 3 prime, keep reads >100bp
allc_path5 = '04_output/trim_r1_min_100bp/bismark_meths/cx_report/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz'
allc5 = pd.read_csv(
    allc_path5,
    delimiter="\t",
    names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
    )
conversion5 = conversion_rate(allc5)

allc5_lambda = allc5.loc[allc5['chr'] == "Lambda_NEB"]
float(allc5_lambda['mC'].sum()) / (allc5_lambda['mC'].sum() + allc5_lambda['uC'].sum())

for c in allc5['chr'].unique():
    chr = allc5.loc[allc5['chr'] == c]
    print(c, chr['mC'].sum(), chr['uC'].sum())



# Trim 5 ad 3 prime by 50 bp
allc_path6 = '04_output/trim_50/bismark_meths/cx_report/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz'
allc6 = pd.read_csv(
    allc_path6,
    delimiter="\t",
    names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
    )
conversion6 = conversion_rate(allc6)

# Specificying libraries should be directional
allc_path7 = '04_output/directional/bismark_meths/cx_report/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz'
allc7 = pd.read_csv(
    allc_path7,
    delimiter="\t",
    names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
    )
conversion7 = conversion_rate(allc7)

allc7.loc[(allc7['chr'] == "Lambda_NEB") & (allc7['mC'] > 0 )]


# Specificying libraries should be directional
allc_path8 = '04_output/no_dovetail/bismark_meths/cx_report/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz'
allc8 = pd.read_csv(
    allc_path8,
    delimiter="\t",
    names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
    )
conversion8 = conversion_rate(allc8)


allc7.loc[allc7['chr'] == "Lambda_NEB"]\
    .loc[ (allc7['mC'] >0) | (allc7['uC'] > 0 )]

allc7.loc[allc7['chr'] == "ChrC"].loc[ (allc7['mC'] >0) ]['uC'].sum()
allc1.loc[allc1['chr'] == "ChrC"]['uC'].sum()

# Col-0 control run with pentuple mutant
allc_path_col0 = '04_output/mutants/bismark_meths/cx_report/211006_TGCCGGTCAGCCTGATACAA_S1_L001_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz'
allc_col0 = pd.read_csv(
    allc_path_col0,
    delimiter="\t",
    names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
    )
conversion_col0 = conversion_rate(allc_col0)


from glob import glob
mutant_files = glob('04_output/mutants/bismark_meths/cx_report/*.CX_report.txt.gz')
import os

mutant_mC = {}
for f in mutant_files:  
    print( os.path.basename(f) )
    allc_mutants = pd.read_csv(
        f,
        delimiter="\t",
        names=['chr', 'pos','strand', "mC", "uC", "context", "seq"]
        )
    allc_mutants = allc_mutants.loc[allc_mutants['chr'] == "Lambda_NEB"]
    mutant_mC[os.path.basename(f)] =  float(allc_mutants['mC'].sum() ) / (allc_mutants['mC'].sum() + allc_mutants['uC'].sum() )

pd.DataFrame(mutant_mC)


# An exmaple from the old protocol
allc_path_old = '/groups/nordborg/projects/epiclines/003.dogging_expt/004.output/001.methylseq/methylpy/allc_HY5W7DMXX_1#148028_TAGCGCTCGCGTATAT.tsv.gz'
allc_old = pd.read_csv(
    allc_path_old,
    delimiter="\t",
    names=['chr', 'pos','strand', "context", "mC", "total", "call"]
    )
conversion_old = allc_old.groupby('chr').apply(
    lambda x: float(x['mC'].sum()) / x['total'].sum()
)
