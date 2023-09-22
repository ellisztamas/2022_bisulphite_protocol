#!/bin/bash

# Call methylation status in 150bp windows across the each genome 

#SBATCH --job-name=mC_windows_col0_flies
#SBATCH --qos=medium
#SBATCH --time=1-00:00:00
#SBATCH --mem=20gb
#SBATCH --output=slurm/mC_windows_col0_flies.%a.out
#SBATCH --error=slurm/mC_windows_col0_flies.%a.err
#SBATCH --array=0-11

module load build-env/f2022
module load anaconda3/2023.03
source ~/.bashrc
conda activate epiclines

date

# Array of four cytosine reports
# files=(03_processing/03_map_original_f2s/output/reports/*CX_report.txt.gz)
files=(/scratch-cbe/users/thomas.ellis/05_col0_flies/methylseq/bismark/coverage2cytosine/reports/*CX_report.txt.gz)

# File for this job
infile=${files[$SLURM_ARRAY_TASK_ID]}
# Output directory
outdir=03_processing/03_map_original_f2s/output/windows/
mkdir -p $outdir
# output filename for methylation in windows
outfile_windows=$outdir/$(basename $infile)
outfile_windows=${outfile_windows/CX_report.txt.gz/mC_in_windows.csv}

python 02_library/python/mC_in_windows.py \
--input $infile \
--output $outfile_windows \
--window 150

date