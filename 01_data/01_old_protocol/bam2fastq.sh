#!/bin/bash

# Tom Ellis
# Convert BAM files to fastq

# SLURM
#SBATCH --mem=5GB
#SBATCH --output=slurm/bam2fastqc.%J.out
#SBATCH --error=slurm/bam2fastqc.%J.err
#SBATCH --qos=rapid
#SBATCH --time=10:00
#SBATCH --array=0-4

ml build-env/f2021
ml bedtools/2.30.0-gcc-10.2.0

indir=/groups/nordborg/projects/epiclines/003.dogging_expt/001.data/003.sequencing/mergeAP
files=($indir/*bam)
infile=${files[$SLURM_ARRAY_TASK_ID]}

outdir=01_data/01_old_protocol
outfile=`basename ${infile/.bam/.fastq}`

echo $infile
echo $outfile

bedtools bamtofastq \
-i  $infile \
-fq $outdir/$outfile