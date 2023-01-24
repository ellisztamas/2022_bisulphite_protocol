#!/bin/bash
#SBATCH --job-name=bam2sam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --qos=rapid
#SBATCH --mem=10gb
#SBATCH --array=0-4
#SBATCH --output=./slurm/meth_counts.%J.out
#SBATCH --error=./slurm/meth_counts.%J.err

module load build-env/f2022
module load samtools/1.15-gcc-11.2.0

# export PYTHONPATH="/groups/nordborg/projects/epiclines/006.quality_control/01_2022_bisulphite_protocol/02_library"

files=(05_results/05_disembowel_samfiles/samfiles/*sam)
infile=${files[$SLURM_ARRAY_TASK_ID]}

outdir=05_results/05_disembowel_samfiles/meth_counts
mkdir -p $outdir
outfile=`basename ${infile/.sam/.tsv}`

python 02_library/python/BismarkSam.py \
-i $infile \
-o $outdir/$outfile