#!/bin/bash

# Script to sort the deduplicated BAM files from the methylseq pipeline
# Tom Ellis, 12th September 2023

#SBATCH --job-name=sort_col0_flies
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --qos=rapid
#SBATCH --mem=10gb
#SBATCH --array=0-12
#SBATCH --output=slurm/sort_col0_flies.%a.out
#SBATCH --error=slurm/sort_col0_flies.%a.err

module load build-env/f2022
module load anaconda3/2023.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate epiclines

scratch=/scratch-cbe/users/$(whoami)/05_col0_flies/methylseq/bismark

files=($scratch/deduplicated/*.bam)
infile=${files[$SLURM_ARRAY_TASK_ID]}

outfile=${infile//deduplicated/sorted}
mkdir -p $(dirname $outfile)

samtools sort -o $outfile $infile
samtools index ${outfile}

