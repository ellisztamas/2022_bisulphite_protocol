# Create a sample sheet giving locations of raw fastq files for the original 
# samples that were subsequently repeated on the mix plate.

# This binds and subsets the (much larger) sample sheets for (most) parents and F2s
# from Rahul Pisupati's experiment. See here for how that was done:
# /groups/nordborg/projects/epiclines/007.hmm/03_processing/03_align_plate_positions/

# Note that not all the original samples are available, because there are two
# plates from the original sample where we are not 100% sure we have the right data

# Tom Ellis,7th September 2023

import pandas as pd

# File giving the experimental design with a ow for each sample, giving plate 
# ID, temperature, row and column, F2 vs parent, and direction.
design = pd.read_csv("01_data/12_resequenced_f2s/experimental_design.csv")
design['design_sample_id'] = design['sample'] + "_" + design['genotype'].astype(str)

# Import the two sample sheets previously created for F2s and parent samples
# Merge into a single dataframe
original_sample_sheet = pd.concat([
    pd.read_csv("03_processing/03_map_original_f2s//f2_plate_positions.csv"),
    pd.read_csv("03_processing/03_map_original_f2s//parents_plate_positions.csv")
])

# Subset the dataframe to return only those samples repeated in the mix plate
subsetted_sample_sheet = original_sample_sheet.loc[original_sample_sheet['sample'].isin(design['design_sample_id'])]

# Save to disk
subsetted_sample_sheet.to_csv(
    "03_processing/03_map_original_f2s/original_plate_positions.csv", index = False
)