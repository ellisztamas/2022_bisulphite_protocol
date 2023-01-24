#!/bin/bash
#SBATCH --job-name=meth_counts
#SBATCH --nodes=1
#SBATCH --qos=rapid
#SBATCH --mem=10gb
#SBATCH --array=0-3
#SBATCH --output=./slurm/meth_counts.%J.out
#SBATCH --error=./slurm/meth_counts.%J.err

# Convert files from BAM to SAM
module load build-env/f2022
module load samtools/1.15-gcc-11.2.0
# Convert BAM files from the spike in with the new protocol
files=(04_output/three_prime_15/sorted_bam/*.bam)
bamfile=${files[$SLURM_ARRAY_TASK_ID]}
# Output
tempdir=05_results/05_disembowel_samfiles/samfiles
mkdir -p $tempdir
samfile=`basename ${bamfile/.bam/.sam}`
# Conversion
samtools view $bamfile > $tempdir/$samfile

# Extract methylation counts
export PYTHONPATH="/groups/nordborg/projects/epiclines/006.quality_control/01_2022_bisulphite_protocol/02_library/python"
# Input 
samfiles=(05_results/05_disembowel_samfiles/samfiles/*sam)
infile=${samfiles[$SLURM_ARRAY_TASK_ID]}
# Output
outdir=05_results/05_disembowel_samfiles/meth_counts
mkdir -p $outdir
outfile=`basename ${infile/.sam/.tsv}`
# Run the script to
python 05_results/05_disembowel_samfiles/meth_counts.py \
-i $infile \
-o $outdir/$outfile

rm $tempdir/$samfile