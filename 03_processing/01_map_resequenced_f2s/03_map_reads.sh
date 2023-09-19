# Code to map bisulphite short reads to the TAIR reference genome.
# This uses the methylseq pipeline (https://nf-co.re/methylseq/2.4.0)

# The pipeline requires raw data to be in fastq.gz format, and a sample sheet stating
# how sample names line up with files. These are created with other scripts in 
# this folder. 

# Tom Ellis, 3rd July 2023

module load build-env/f2022
module load nextflow/22.10.7

# set the working directory on scratch
scratch=/scratch-cbe/users/$(whoami)/mix_plate
mkdir -p $scratch

# Location of the methylseq nextflow pipeline
pipeline_location=02_library/methylseq/main.nf
# Map reads to a Fasta file with the two parents plus TAIR10 organelles
ref_genome="01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta"

nextflow -log $scratch/nextflow.log run $pipeline_location \
--input 03_processing/01_map_resequenced_f2s/mix_plate_positions.csv \
--outdir $scratch/methylseq \
--fasta $ref_genome \
--clip_r1 15 --clip_r2 15  --three_prime_clip_r1 9  --three_prime_clip_r2 9 \
--non_directional \
--cytosine_report \
--aligner bismark \
-profile cbe \
-w $scratch/work

