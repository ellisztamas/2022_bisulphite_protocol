#!/bin/bash

# Script to sort the deduplicated BAM files from the methylseq pipeline
# Tom Ellis, 12th September 2023

#SBATCH --job-name=sort_originals
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --qos=rapid
#SBATCH --mem=10gb
#SBATCH --array=0-59
#SBATCH --output=slurm/sort_originals.%a.out
#SBATCH --error=slurm/sort_originals.%a.err

module load build-env/f2022
module load samtools/1.15-gcc-11.2.0

scratch=/scratch-cbe/users/$(whoami)/original_f2s/methylseq/bismark

files=($scratch/deduplicated/*.bam)
infile=${files[$SLURM_ARRAY_TASK_ID]}

outfile=${infile//deduplicated/sorted}
mkdir -p $(dirname $outfile)

samtools sort -o $outfile $infile
samtools index ${outfile}

