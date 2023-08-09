#!/bin/bash

# Script to extract only reads mapping to the chloroplast for a single file

#SBATCH --job-name=extract_chloroplast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --qos=rapid
#SBATCH --mem=10gb
#SBATCH --output=./slurm/extract_chloroplast.%a.out
#SBATCH --error=./slurm/extract_chloroplast.%a.err

module load build-env/f2022
module load samtools/1.15-gcc-11.2.0

infile=04_output/three_prime_15/sorted_bam/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam
outdir=05_results/05_disembowel_samfiles/samfiles
mkdir -p $outdir

samtools view $infile "ChrC" > $outdir/chloroplast.sam