#!/bin/bash
#SBATCH --job-name=discard_reads
#SBATCH --nodes=1
#SBATCH --qos=rapid
#SBATCH --time=30:00
#SBATCH --mem=50gb
#SBATCH --array=0-4
#SBATCH --output=./slurm/discard_reads.%J.out
#SBATCH --error=./slurm/discard_reads.%J.err

module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate bisulphite-qc

# Working subdirectory
dir=05_results/06_discard_clustered_reads

# Convert files from BAM to SAM
# Convert BAM files from the spike in with the new protocol
files=(04_output/three_prime_15/bismark_dedup/*.bam)
bamfile=${files[$SLURM_ARRAY_TASK_ID]}
# Output
mkdir -p $dir/samfiles
samfile=`basename ${bamfile/.deduplicated.bam/.sam}`
# BAM to SAM Conversion
samtools view -f 0x2 -h $bamfile > $dir/samfiles/$samfile

# Split reads based on whether they have cluster of cytosines or not
export PYTHONPATH="/groups/nordborg/projects/epiclines/006.quality_control/01_2022_bisulphite_protocol/02_library/python"
# Output
mkdir -p $dir/cleaned_samfiles
mkdir -p $dir/filtered_reads
# Run the script
python $dir/discard_clustered_reads.py \
--input $dir/samfiles/$samfile \
--keep $dir/cleaned_samfiles/$samfile \
--discard $dir/filtered_reads/$samfile

# Extract methylation on each read
bismark_methylation_extractor \
$dir/cleaned_samfiles/$samfile \
--genome_folder 01_data/03_reference_genome \
--paired-end \
--no_overlap \
--comprehensive \
--bedGraph \
--cytosine_report \
--mbias_off \
--output_dir $dir/bismark_meths

# Get a single cytosine report
covfile=`basename ${samfile/.sam/.bismark.cov.gz}`
outfile=`basename ${samfile/.sam/}`
coverage2cytosine \
$dir/bismark_meths/$covfile \
--gzip \
--genome 01_data/03_reference_genome \
--CX_context \
--dir 05_results/06_discard_clustered_reads/bismark_meths/ \
-o $outfile

# Create and index sorted BAM files you can view in IGV
mkdir -p $dir/sorted_bams
sortedfile=`basename ${bamfile/.deduplicated.bam/.sorted.bam}`
samtools view -S -h -b $dir/cleaned_samfiles/$samfile | samtools sort -m 100M -o $dir/sorted_bams/$sortedfile
samtools index -b $dir/sorted_bams/$sortedfile -o $dir/sorted_bams/$sortedfile.bai

mkdir -p $dir/sorted_bams_discard
samtools view -S -h -b $dir/filtered_reads/$samfile | samtools sort -m 100M -o $dir/sorted_bams_discard/$sortedfile
samtools index -b $dir/sorted_bams_discard/$sortedfile -o $dir/sorted_bams_discard/$sortedfile.bai


# Remove the intermediate SAM file
rm $dir/samfiles/$samfile
rm $dir/cleaned_samfiles/$samfile
rm $dir/filtered_reads/$samfile
rm $dir/bismark_meths/*txt
