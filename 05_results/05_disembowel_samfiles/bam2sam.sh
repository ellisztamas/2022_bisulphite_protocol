#!/bin/bash
#SBATCH --job-name=bam2sam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --qos=rapid
#SBATCH --mem=10gb
#SBATCH --array=0-4
#SBATCH --output=./slurm/bam2sam.%J.out
#SBATCH --error=./slurm/bam2sam.%J.err

module load build-env/f2022
module load samtools/1.15-gcc-11.2.0

files=(04_output/three_prime_15/sorted_bam/*.bam)
infile=${files[$SLURM_ARRAY_TASK_ID]}

outdir=05_results/05_disembowel_samfiles/samfiles
mkdir -p $outdir
outfile=`basename ${infile/.bam/.sam}`

samtools view $infile > $outdir/$outfile