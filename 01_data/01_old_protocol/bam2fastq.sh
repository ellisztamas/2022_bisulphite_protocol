#!/bin/bash

# Tom Ellis
# Convert BAM files to fastq

# SLURM
#SBATCH --mem=5GB
#SBATCH --output=slurm/bam2fastqc.%J.out
#SBATCH --error=slurm/bam2fastqc.%J.err
#SBATCH --qos=medium
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

ml build-env/f2021
ml bedtools/2.30.0-gcc-10.2.0
# ml build-env/f2022
# ml samtools/1.15-gcc-11.2.0

indir=/resources/ngs/nordborg/12091
outdir=01_data/01_old_protocol

# FILES=($indir/*bam)
# infile=$FILES[${SLURM_ARRAY_TASK_ID}]

# filename=`basename $infile`
# # outfile=$outdir/${filename%.*}.fastq

# echo $infile
# echo $filename
# echo $outfile

bedtools bamtofastq \
-i   $indir/HHHY7DSX2_2#169198_ACTCGCTAAAGGAGTA.bam \
-fq $outdir/HHHY7DSX2_2#169198_ACTCGCTAAAGGAGTA.fastq

bedtools bamtofastq \
-i   $indir/HHHY7DSX2_2#169198_ACTCGCTAACTGCATA.bam \
-fq $outdir/HHHY7DSX2_2#169198_ACTCGCTAACTGCATA.fastq

bedtools bamtofastq \
-i   $indir/HHHY7DSX2_2#169198_ACTCGCTAAGAGTAGA.bam \
-fq $outdir/HHHY7DSX2_2#169198_ACTCGCTAAGAGTAGA.fastq