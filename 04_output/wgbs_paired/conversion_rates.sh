#!/bin/bash
#SBATCH --job-name=conversion_rates
#SBATCH --nodes=1
#SBATCH --qos=rapid
#SBATCH --time=10:00
#SBATCH --mem=10gb
#SBATCH --array=0-4
#SBATCH --output=./slurm/conversion_rates.%J.out
#SBATCH --error=./slurm/conversion_rates.%J.err

module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate bisulphite-qc

files=(04_output/wgbs_paired/bismark_meths/cx_report/*)
infile=${files[$SLURM_ARRAY_TASK_ID]}

outfile=04_output/wgbs_paired/conversion_rates.csv

python 02_library/python/calculate_conversion.py \
--allc $infile \
--output $outfile
