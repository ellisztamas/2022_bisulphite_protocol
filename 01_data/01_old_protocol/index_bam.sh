#!/bin/sh

# Tom Ellis
# Sort Bam file and create and index file.

# SLURM
#SBATCH --mem=5GB
#SBATCH --output=./01_data/04_pentuple_mutant/log
#SBATCH --qos=medium
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

ml build-env/f2022
ml samtools/1.15-gcc-11.2.0

samtools sort 01_data/01_old_protocol/mix_P6.A10.bam -o 01_data/01_old_protocol/mix_P6.A10_sorted.bam
samtools index -b 01_data/01_old_protocol/mix_P6.A10_sorted.bam > 01_data/01_old_protocol/mix_P6.A10_sorted.bai
