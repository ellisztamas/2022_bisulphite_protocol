#!/bin/bash
#SBATCH --job-name=fastqc_lambdaNEB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --mem=10gb
#SBATCH --output=slurm/fastqc_lambdaNEB.%a.out
#SBATCH --error=slurm/fastqc_lambdaNEB.%a.err

input_bam='05_results/02_lambda_reads/output/filtered_bams/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam'

ml build-env/f2021
ml fastqc/0.11.9-java-11

outdir=05_results/02_lambda_reads/output/fastqc
mkdir -p $outdir

fastqc $input_bam -o $outdir