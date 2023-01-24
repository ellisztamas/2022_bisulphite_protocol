#!/bin/bash
#SBATCH --job-name=SAMtools_View
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --mem=10gb
#SBATCH --output=./slurm/SAMtools.%J.out
#SBATCH --error=./slurm/SAMtools.%J.err

module load build-env/f2022
module load samtools/1.15-gcc-11.2.0

input_bam='04_output/wgbs_paired/sorted_bam/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam'

outdir=05_results/02_lambda_reads/output/filtered_bams
mkdir -p $outdir

# samtools view -bhS -@ $SLURM_NTASKS_PER_NODE $input_bam -o $outdir/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam

samtools view -h $input_bam "Lambda_NEB:3000-3300" > $outdir/test.bam