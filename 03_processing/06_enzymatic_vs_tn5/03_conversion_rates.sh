#!/bin/bash

# Calculate methylation levels on each chromosome.
# Tom Ellis, 22nd September 2023

#SBATCH --job-name=conversion_em_tn5
#SBATCH --time=1:00:00
#SBATCH --qos=rapid
#SBATCH --mem=10gb
#SBATCH --array=0-11 # Start at 1, because the sample sheet has a header row
#SBATCH --output=slurm/conversion_em_tn5-%a.out
#SBATCH --error=slurm/conversion_em_tn5-%a.err

module load build-env/f2022
module load anaconda3/2023.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate epiclines

indir=03_processing/06_enzymatic_vs_tn5/output

files=($indir/reports/*CX_report.txt.gz)
infile=${files[$SLURM_ARRAY_TASK_ID]}

outdir=$indir/conversion
mkdir -p $outdir

outfile=$(basename ${infile/CX_report.txt.gz/conversion.csv})

python 02_library/python/conversion_rate.py \
    --input $infile \
    --output $outdir/$outfile