#!/bin/bash

# Calculate methylation levels on each chromosome.
# Tom Ellis, 22nd September 2023

#SBATCH --job-name=conversion_mixplate
#SBATCH --output=slurm/conversion_mixplate-%a.out
#SBATCH --error=slurm/conversion_mixplate-%a.err
#SBATCH --time=1:00:00
#SBATCH --qos=rapid
#SBATCH --mem=10gb
#SBATCH --array=0-95

module load build-env/f2022
module load anaconda3/2023.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate epiclines

indir=03_processing/02_map_resequenced_f2s/output

files=($indir/reports/*CX_report.txt.gz)
infile=${files[$SLURM_ARRAY_TASK_ID]}

outdir=$indir/conversion
mkdir -p $outdir

outfile=$(basename ${infile/CX_report.txt.gz/conversion.csv})

python 02_library/python/conversion_rate.py \
    --input $infile \
    --output $outdir/$outfile
