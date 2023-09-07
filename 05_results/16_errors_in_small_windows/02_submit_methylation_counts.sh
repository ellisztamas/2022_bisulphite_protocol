#!/bin/bash

# Submit the python script 01_get_methylation_counts.py as a SLURM job to get 
# non-conversion rate for each of four cytosine report files.
# Tom Ellis

#SBATCH --job-name=16_small_windows
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --mem=20gb
#SBATCH --output=slurm/16_small_windows.%a.out
#SBATCH --error=slurm/16_small_windows.%a.err
#SBATCH --array=0-3

module load build-env/f2022
module load anaconda3/2023.03
source ~/.bashrc
conda activate epiclines

date

# Array of four cytosine reports
files=(04_output/30x_col0/bismark_meths/cx_report/*CX_report.txt.gz)

# File for this job
infile=${files[$SLURM_ARRAY_TASK_ID]}
# Output directory
outdir=05_results/16_errors_in_small_windows/output
mkdir -p $outdir
# output filename for methylation in windows
outfile_windows=$outdir/$(basename $infile)
outfile_windows=${outfile_windows/CX_report.txt.gz/mC_in_windows.csv}

python 05_results/16_errors_in_small_windows/01_get_methylation_counts.py \
--input $infile \
--windows $outfile_windows
date