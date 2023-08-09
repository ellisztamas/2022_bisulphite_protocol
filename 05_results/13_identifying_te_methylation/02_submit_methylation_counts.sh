#!/bin/bash

# Submit the python script 01_get_methylation_counts.py as a SLURM job
# Tom Ellis

#SBATCH --job-name=11_get_methylation_counts
#SBATCH --qos=short
#SBATCH --time=8:00:00
#SBATCH --mem=80gb
#SBATCH --output=slurm/11_get_methylation_counts.out
#SBATCH --error=slurm/11_get_methylation_counts.err

module load build-env/f2022
module load anaconda3/2023.03
source ~/.bashrc
conda activate epiclines

date

python 05_results/13_identifying_te_methylation/01_get_methylation_counts.py

date