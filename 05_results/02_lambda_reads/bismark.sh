#!/bin/bash
#SBATCH --job-name=bismark_lambda
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --mem=10gb
#SBATCH --output=./slurm/bismark_lambda.%a.out
#SBATCH --error=./slurm/bismark_lambda.%a.err

input_fq='05_results/02_lambda_reads/output/filtered_bams/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001_val_1_bismark_bt2_pe.sorted.fastq'
genome=01_data/03_reference_genome
outdir=05_results/02_lambda_reads/output/bismark
mkdir -p $outdir

module load build-env/2020
module load bismark/0.22.2-foss-2018b

bismark \
--parallel 8 \
--genome_folder $genome \
--non_directional \
-o $outdir \
$input_fq