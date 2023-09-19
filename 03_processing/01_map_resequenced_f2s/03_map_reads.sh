#!/bin/bash

# Map bisulphite short reads to the TAIR reference genome.
#
# The pipeline requires raw data to be in fastq.gz format, and a sample sheet stating
# how sample names line up with files. These are created with other scripts in 
# this folder. 
#
# Note that the script hard codes various optional arguments to trim galore and
# Bismark - check that script for which arguments are used.

# Tom Ellis, 19th September 2023

#SBATCH --job-name=align_mixplate
#SBATCH --time=2-00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5G
#SBATCH --qos=medium
#SBATCH --array=1-96 # Start at 1, because the sample sheet has a header row
#SBATCH --output=slurm/align_mixplate-%a.out
#SBATCH --error=slurm/align_mixplate-%a.err

module load build-env/f2022
module load anaconda3/2023.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate epiclines

# working directory
scratch=/scratch-cbe/users/thomas.ellis/mix_plate
# output directory
outdir=03_processing/01_map_resequenced_f2s/output
# FASTA file for the genome to map to
genome=01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta

# Sample sheet giving sample name and paths to the two fastq files
sample_sheet=03_processing/01_map_resequenced_f2s/mix_plate_positions.csv
# Get the sample name
sample_names=$(cut -d',' -f1 $sample_sheet)
sample_names=($sample_names)
# Path for read 1
read1_col=$(cut -d',' -f2 $sample_sheet)
read1_col=($read1_col)
# Path for read 2
read2_col=$(cut -d',' -f3 $sample_sheet)
read2_col=($read2_col)

# Run the script
02_library/bash/bismark_pipeline.sh \
    --sample "${sample_names[$SLURM_ARRAY_TASK_ID]}" \
    --read1  "${read1_col[$SLURM_ARRAY_TASK_ID]}" \
    --read2  "${read2_col[$SLURM_ARRAY_TASK_ID]}" \
    --genome $genome \
    --work $scratch \
    --outdir $outdir