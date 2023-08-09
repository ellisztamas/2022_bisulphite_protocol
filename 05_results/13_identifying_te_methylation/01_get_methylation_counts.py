"""
Script to import the TAIR10 annotation file and pull out information methylated
and unmethylated reads in each context for:
1. genes
2. RdDM- and CMT-regulated TEs
3. windows of 1000 base pairs across each chromosome
for the genotyping file on the 30x Col-0 data we ran.
"""

import pandas as pd
import epiclinestools as epi

print("Using epiclinestools verion {}".format(epi.__version__))

# Cytosine coverage file
print("Importing coverage file.")
path = "04_output/30x_col0/bismark_meths/cx_report/220842_ATGTTGTTGGCAATCTATGA_S9_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz"
col0 = epi.CytosineCoverageFile(path)


# METHYLATION ON GENE BODIES
# annotated genes in tair10
print("Count methylated reads over gene bodies.")
tair10_genes = pd.read_csv(
    "01_data/11_tair10/TAIR10_GFF3_genes.gff", sep = "\t",
    names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    )
tair10_genes = tair10_genes.loc[tair10_genes['type'] == "gene"]
# pull gene names and chromosome labels from the attributes column
tair10_genes['seqid'] = tair10_genes['attributes'].str.slice(3,12)
tair10_genes['chr']   = "Chr" + tair10_genes['seqid'].str.slice(2,3)

gbm = col0.methylation_over_features(
    chr = tair10_genes['chr'],
    start = tair10_genes['start'],
    stop = tair10_genes['end'],
    names = tair10_genes['seqid']
)

gbm.to_csv("05_results/13_identifying_te_methylation/output/mC_genes.csv", index=False)

# CMT2- AND RdDM-REGULATED TRANSPOSABLE ELEMENTS
# Full list of TEs from the TAIR10 annotation
print("Counting methylated reads over transposable elements.")
tair10_TEs = pd.read_csv("01_data/11_tair10/TAIR10_Transposable_Elements.txt", sep = "\t")
tair10_TEs['chr'] = "Chr" + tair10_TEs['Transposon_Name'].str.slice(2,3)

tair10_TEs['Transposon_max_End'] - tair10_TEs['Transposon_min_Start']


# Merge with files for TEs regulated by CMT2 and RdDM
CMT2_TEs = tair10_TEs.merge(
    pd.read_csv("01_data/11_tair10/CMT2_target_TEs.txt", names=['Transposon_Name'])
)
RdDM_TEs = tair10_TEs.merge(
    pd.read_csv("01_data/11_tair10/RdDM_target_TEs.txt", names=['Transposon_Name'])
)
# Get methylation on each TE.
cmt2_meth = col0.methylation_over_features(
    chr = CMT2_TEs['chr'],
    start = CMT2_TEs['Transposon_min_Start'],
    stop = CMT2_TEs['Transposon_max_End'],
    names = CMT2_TEs['Transposon_Name']
)
rddm_meth = col0.methylation_over_features(
    chr = RdDM_TEs['chr'],
    start = RdDM_TEs['Transposon_min_Start'],
    stop = RdDM_TEs['Transposon_max_End'],
    names = RdDM_TEs['Transposon_Name']
)
# Write to disk
cmt2_meth.to_csv("05_results/13_identifying_te_methylation/output/mC_on_cmt2_TEs.csv", index = False)
rddm_meth.to_csv("05_results/13_identifying_te_methylation/output/mC_on_rddm_TEs.csv", index = False)

# WINDOWS
# Get counts in windows of 1000 bp
print("Counting methylated reads in windows")
col0.methylation_in_windows(window_size = 1000).to_csv(
    "05_results/13_identifying_te_methylation/output/mC_in_windows.csv", index=False
)
