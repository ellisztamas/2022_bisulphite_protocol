#!/usr/bin/env bash
# 
# Tom Ellis
# Script to unzip each raw data file to the `scratch-cbe` drive
# on the VBC cluster.

# SLURM
#SBATCH --mem=5GB
#SBATCH --output=./03_processing/unzip_raw_bams.log
#SBATCH --qos=medium
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8

# Reference directories. Change these for your machine.
# Home folder for the project
PROJ=/groups/nordborg/projects/epiclines/006.quality_control/01_2022_bisulphite_protocol
# Working directory to perform computations on the VBC cluster.
DIR=/scratch-cbe/users/thomas.ellis

# Where the data are
# Location of the raw zip file on the VBC cluster
RAW=$PROJ/01_data/02_new_protocol/BRL95611-0924_0_R13786_20220701.tar.gz

# Where to save the output
DATA=$DIR/01_unzipped_raw_bams # where to unzip raw reads
mkdir -p $DATA

# Unzip raw data
tar -xvC $DATA -f $RAW --wildcards '*.fastq.gz'