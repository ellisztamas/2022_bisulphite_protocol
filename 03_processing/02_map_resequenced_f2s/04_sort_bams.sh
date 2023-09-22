#!/bin/bash

# Script to sort the deduplicated BAM files from the methylseq pipeline
# Tom Ellis, 12th September 2023

#SBATCH --job-name=sort_mix_plate
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --qos=rapid
#SBATCH --mem=10gb
#SBATCH --array=0-95
#SBATCH --output=slurm/sort_mix_plate.%a.out
#SBATCH --error=slurm/sort_mix_plate.%a.err

module load build-env/f2022
module load samtools/1.15-gcc-11.2.0

scratch=/scratch-cbe/users/$(whoami)/mix_plate/methylseq/bismark

files=($scratch/deduplicated/*.bam)
infile=${files[$SLURM_ARRAY_TASK_ID]}

outfile=${infile//deduplicated/sorted}
mkdir -p $(dirname $outfile)

samtools sort -o $outfile $infile
samtools index ${outfile}

