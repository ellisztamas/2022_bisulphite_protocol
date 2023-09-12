#!/bin/bash

# Call methylation status in 150bp windows across the each genome 

#SBATCH --job-name=mC_windows_originals
#SBATCH --qos=medium
#SBATCH --time=1-00:00:00
#SBATCH --mem=20gb
#SBATCH --output=slurm/mC_windows_originals.%a.out
#SBATCH --error=slurm/mC_windows_originals.%a.err
#SBATCH --array=0-58

module load build-env/f2022
module load anaconda3/2023.03
source ~/.bashrc
conda activate epiclines

date

# Array of four cytosine reports
files=(03_processing/02_map_original_f2s/output/reports/*CX_report.txt.gz)

# File for this job
infile=${files[$SLURM_ARRAY_TASK_ID]}
# Output directory
outdir=03_processing/02_map_original_f2s/output/windows/
mkdir -p $outdir
# output filename for methylation in windows
outfile_windows=$outdir/$(basename $infile)
outfile_windows=${outfile_windows/CX_report.txt.gz/mC_in_windows.csv}

python 02_library/python/mC_in_windows.py \
--input $infile \
--output $outfile_windows \
--window 1000

date