#!/bin/bash
#SBATCH --job-name=which_reads_fail
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --mem=20gb
#SBATCH --output=./slurm/which_read_fails.py.%a.out
#SBATCH --error=./slurm/which_read_fails.%J.err

conda activate pybshap

python 05_results/01_which_reads_fail/which_reads_fail.py