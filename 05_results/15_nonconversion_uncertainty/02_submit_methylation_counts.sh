#!/bin/bash

# Submit the python script 01_get_methylation_counts.py as a SLURM job to get 
# non-conversion rate for each of four cytosine report files.
# Tom Ellis

#SBATCH --job-name=15_nonconversion
#SBATCH --qos=short
#SBATCH --time=8:00:00
#SBATCH --mem=80gb
#SBATCH --output=slurm/15_nonconversion.%a.out
#SBATCH --error=slurm/15_nonconversion.%a.err
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
outdir=05_results/15_nonconversion_uncertainty/output/col0
mkdir -p $outdir
# output filename for methylation in windows
outfile_windows=$outdir/$(basename $infile)
outfile_windows=${outfile_windows/CX_report.txt.gz/mC_in_windows.csv}
# output filename for methylation over annotated features
outfile_features=$outdir/$(basename $infile)
outfile_features=${outfile_features/CX_report.txt.gz/mC_over_features.csv}

python 05_results/15_nonconversion_uncertainty/01_get_methylation_counts.py \
--input $infile \
--windows $outfile_windows \
--features $outfile_features

date