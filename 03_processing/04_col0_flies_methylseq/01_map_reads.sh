# Code to map bisulphite short reads for Columbia and Drosophila samples sampled
# at high coverage using the methylseq pipeline (https://nf-co.re/methylseq/2.4.0)
# Raw data are unzipped and renamed in 01_data/05_col0_flies/unzip_tarball.sh

# The pipeline requires raw data to be in fastq.gz format, and a sample sheet stating
# how sample names line up with files.

# Tom Ellis, 14th September 2023

module load build-env/f2022
module load nextflow/22.10.7

# set the working directory on scratch
scratch=/scratch-cbe/users/$(whoami)/05_col0_flies
mkdir -p $scratch
# Location of the methylseq nextflow pipeline
pipeline_location=02_library/methylseq/main.nf

# # Columbia samples
# # Map reads to a the TAIR10 genome including organelles and the Lambda and pUC19 vectors
# tair10_genome="01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta"
# nextflow -log $scratch/col0/nextflow.log run $pipeline_location \
# --input 03_processing/04_col0_flies_methylseq/col0_sample_sheet.csv \
# --outdir $scratch/methylseq \
# --fasta $tair10_genome \
# --clip_r1 15 --clip_r2 15  --three_prime_clip_r1 9  --three_prime_clip_r2 9 \
# --non_directional \
# --cytosine_report \
# --aligner bismark \
# -profile cbe \
# -w $scratch/col0/work

# Flies
# Map to the drosophila 6 genome, including bacterial genomes and Lambda and pUC19 vectors
drosophila_genome=01_data/06_fly_genome/dm6_plus_vectors_microbiome.fasta
nextflow -log $scratch/d_mel/nextflow.log run $pipeline_location \
--input 03_processing/04_col0_flies_methylseq/drosophila_sample_sheet.csv \
--outdir $scratch/methylseq \
--fasta $drosophila_genome \
--clip_r1 15 --clip_r2 15  --three_prime_clip_r1 9  --three_prime_clip_r2 9 \
--non_directional \
--cytosine_report \
--aligner bismark \
-profile cbe \
-w $scratch/d_mel/work
