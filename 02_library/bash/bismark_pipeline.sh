#!/bin/bash

# Script to trim and map bisulphite reads using Trim Galore! and Bismark, then
# report methylation state across the genome using bismark_methylation_extractor.
#
# Trim Galore! and Bismark are run with a minimum of arguments. Additional options
# can (and probably should) be included with the `trim_galore_args` and 
# `bismark_args` options. Be very careful that these options are surrounded by "".

# Tom Ellis 19th September 2023

# INPUTS:
# sample: Name of the sample
# read1: Path to a file containing data for read 1 if data are paired end, or 
#     else all read if data are single end.
# read2: Path to a file containing data for read 2 if data are paired end.
# genome: Path to a FASTA genome to map reads to.
#    It is assumed you have previously run bismark_genome_preparation on this
    # genome.
# work: Path to a directory to save files to as the pipeline runs.
# outdir: Path to a directory to copy the output files to. Not files are staged
    # out automatically.
# trim_galore_args: Additional arguments to be passed to Trim Galore! This should
#   be surrounded by "".
# bismark_args: Additional arguments to be passed to Bismark. This should be
#   surrounded by "".

# OUTPUT:
# The script will copy the most useful files to the directory specified in `outdir`:
# outdir
#     - logs
#         - bismark_report.txt   
#         - bismark_deduplication_report.txt
#     - reports
#         CX reports for each samples
#     - sorted
#         Aligned BAM files sorted by position
#         Index files for those BAM files

# EXAMPLE USAGE
# sample_name=mix_A1
# read1=$scratch/fastq/H7MWJDSX3_1#188467_TGCCGGTCAGCCTGATACAA.fastq.gz
# read2=$scratch/fastq/H7MWJDSX3_2#188467_TGCCGGTCAGCCTGATACAA.fastq.gz
# genome=01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fa
# work=/scratch-cbe/users/thomas.ellis/mix_plate/

# 02_library/bash/bismark_pipeline.sh \
# --sample $sample_name \
# --read1 $read1 \
# --read2 $read2 \
# --genome $genome \
# --work $scratch \
# --outdir 03_processing/01_map_resequenced_f2s/output \
# --trim_galore_args "--clip_r1 15 --clip_r2 15 --three_prime_clip_R1 9 --three_prime_clip_R2 9 --cores 4" \
# --bismark_args "--local --non_directional --strandID"

date
echo ""

# By default set no additional options to trim galore or Bismark
trim_galore_args=''
bismark_args=''

# *** Make sure you have a new enough getopt to handle long options (see the man page)
getopt -T &>/dev/null
if [[ $? -ne 4 ]]; then echo "Getopt is too old!" >&2 ; exit 1 ; fi

# Parse command line options using getopt
long_options="help,sample:,read1:,read2:,genome:,work:,outdir:,trim_galore_args:,bismark_args:"
options=$(getopt -o h --long "$long_options" -n "$(basename "$0")" -- "$@")
eval set -- "$options"

# Process command line options
while true; do
  case $1 in
    -h | --help)
      usage
      exit 0
      ;;
    --sample) sample_name=$2 ; shift 2;;
    --read1) read1=$2        ; shift 2;;
    --read2) read2=$2        ; shift 2;;
    --genome) genome=$2      ; shift 2;;
    --work) work=$2          ; shift 2;;
    --outdir) outdir=$2      ; shift 2;;
    --trim_galore_args)
      trim_galore_args="$2"
      shift 2
      ;;
    --bismark_args)
      bismark_args="$2"
      shift 2
      ;;
    --)
      shift
      break
      ;;
    *)
      usage
      exit 1
      ;;
  esac
done


# Print information about the input parameters.
echo "Mapping reads for sample $sample_name"
echo "Path to file with mate 1 reads: $read1"
echo "Path to file with mate 2 reads: $read2"

echo "Mapping reads to: $genome"
echo "Setting working directory to $work"
echo "Using trim_galore arguments: $trim_galore_args"
echo "Using Bismark arguments: $bismark_args"

# Trimming with trim galore!
trim_galore_dir=$work/trimgalore
echo ""
echo "Trimming reads, and saving output to: $trim_galore_dir"
mkdir -p $trim_galore_dir

# trim_galore_cmd="trim_galore \
#     --gzip \
#     --paired \
#     -o ${trim_galore_dir} \
#     --basename ${sample_name} \
#     ${trim_galore_args} \
#     ${read1} ${read2}"
# echo "Command passed to trim galore: $trim_galore_cmd"
# $trim_galore_cmd
# echo "Finished trimming reads."
# date
# echo ""

# echo "Running fastqc on the trimmed reads, and saving to $work/fastqc"
# mkdir -p $work/fastqc
# mkdir -p $work/fastqc/temp
# fastqc \
#   --outdir $work/fastqc
#   --dir $work/fastqc/temp
#   ${trim_galore_dir}/${sample_name}_val_1.fq.gz \
#   ${trim_galore_dir}/${sample_name}_val_2.fq.gz
# echo "Finished quality checks."
# date
# echo ""

# # Align reads to the genome
bismark_dir=$work/bismark/$sample_name
# mkdir -p $bismark_dir
# mkdir -p $bismark_dir/temp
# echo "Aligning reads using Bismark and saving the output to"
# echo "$bismark_dir"

# bismark_cmd="bismark $(dirname ${genome}) \
#     -1 ${trim_galore_dir}/${sample_name}_val_1.fq.gz \
#     -2 ${trim_galore_dir}/${sample_name}_val_2.fq.gz \
#     ${bismark_args} \
#     --temp_dir ${bismark_dir}/temp \
#     --output_dir ${bismark_dir}"
# echo "Command passed to Bismark: ${bismark_cmd}"
# $bismark_cmd
# echo "Finished aligning reads."
# date
# echo ""


# # sort the newly aligned BAM file.
# sorted_by_name=$bismark_dir/$sample_name.sortedByName.bam
# mkdir -p $(dirname $sorted_by_name)
# echo "Sorting the newly aligned BAM file by read name, and saving as:"
# echo $sorted_by_name

# samtools sort \
#     -n \
#     -@ $(nproc) \
#     -o $sorted_by_name \
#     $bismark_dir/$sample_name'_val_1_bismark_bt2_pe.bam' # Need to fill this in!
# echo "Finished sorting."
# date
# echo ""



# # Deduplicate BAM file
# mkdir -p $(dirname $sorted_by_name)
# echo Deduplicating the sorted BAM file and saving as $sample_name.sortedByName.deduplicated.bam

# deduplicate_bismark \
#     --bam \
#     --paired \
#     --outfile $sample_name \
#     -output_dir $bismark_dir \
#     $sorted_by_name

# echo "Finished deduplication."
# date
# echo ""



# # Sort and index the deduplicated file by position
# sorted_by_pos=$bismark_dir/$sample_name.sortedByPos.bam
# mkdir -p $(dirname $sorted_by_pos)
# echo "Sorting the deduplicated BAM file by position, and saving as:"
# echo $sorted_by_pos
# samtools sort \
#     -@ $(nproc) \
#     -o $sorted_by_pos \
#     $bismark_dir/$sample_name.deduplicated.bam
# echo "Finished sorting. Indexing..."
# samtools index  $sorted_by_pos
# echo "Finished indexing."
# date
# echo ""

# extract methylation calls
echo "Running bismark_methylation_extractor, and saving to:"
echo $bismark_dir/reports
genome_folder=$(dirname $genome)

bismark_methylation_extractor \
    --paired-end \
    --cytosine_report \
    --CX_context \
    --genome_folder $(realpath $genome_folder) \
    --no_header \
    --no_overlap \
    --comprehensive \
    --gzip \
    --output_dir $bismark_dir/reports \
    $bismark_dir/$sample_name.deduplicated.bam
echo "Finished bismark_methylation_extractor."
date
echo ""

# Stage out the files to keep
echo "Copying files to the output directory."
mkdir -p $outdir/logs
cp $bismark_dir/*.deduplication_report.txt $outdir/logs
cp $bismark_dir/*_val_1_bismark_bt2_PE_report.txt $output_dir/logs

mkdir -p $outdir/reports
cp $bismark_dir/reports/*.CX_report.txt.gz $outdir/reports
mv $outdir/reports/$sample_name.deduplicated.CX_report.txt.gz $outdir/reports/$sample_name.CX_report.txt.gz

mkdir -p $outdir/sorted
cp $bismark_dir/$sample_name.sortedByPos.bam $outdir/sorted
cp $bismark_dir/$sample_name.sortedByPos.bam.bai $outdir/sorted
echo "Finished copying."
date