#!/bin/bash
#SBATCH --job-name=bamtofastq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --mem=10gb
#SBATCH --output=slurm/bamtofastq.%a.out
#SBATCH --error=slurm/bamtofastq.%a.err

module load build-env/2020
module load bedtools/2.27.1-foss-2018b

input_bam='05_results/02_lambda_reads/output/filtered_bams/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam'

outdir=05_results/02_lambda_reads/output/filtered_bams
mkdir -p $outdir

bedtools bamtofastq -i $input_bam -fq $outdir/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.sorted.fastq

