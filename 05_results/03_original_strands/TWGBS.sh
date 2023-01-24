#!/bin/bash
#SBATCH --job-name=original_strand
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mem=10gb
#SBATCH --output=./slurm/original_strand.py.%J.out
#SBATCH --error=./slurm/original_strand.%J.err

module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate pybshap

indir=/scratch-cbe/users/thomas.ellis/00_temp
outdir=05_results/03_original_strands

python 05_results/03_original_strands/TWGBS_read_pair_reconstruction.py \
--R1_in $indir/200183_AACCAGCCACGCCACAGCAC_S60_L001_R1_001.fastq.gz \
--R2_in $indir/200183_AACCAGCCACGCCACAGCAC_S60_L001_R2_001.fastq.gz \
--R1_out $outdir/read1.fastqz.gz \
--R2_out $outdir/read2.fastqz.gz \
--R1_unassigned $outdir/read1_unassigned \
--R2_unassigned $outdir/read2_unassigned \
--log $outdir/log.txt