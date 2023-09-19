#!/bin/bash

# Script to trim and map bisulphite reads using Trim Galore! and Bismark, then
# report methylation state across the genome using bismark_methylation_extractor.

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
# genome=01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta
# work=/scratch-cbe/users/thomas.ellis/mix_plate/

# 02_library/bash/bismark_pipeline.sh \
# --sample $sample_name \
# --read1 $read1 \
# --read2 $read2 \
# --genome $genome \
# --work $scratch \
# --outdir 03_processing/01_map_resequenced_f2s/output

date
echo ""

# Options passed to trim galore and bismark. Change as necessary
trim_galore_args="--clip_r1 15 --clip_r2 15 --three_prime_clip_R1 9 --three_prime_clip_R2 9 --cores 4"
bismark_args="--local --non_directional --strandID"

# *** Make sure you have a new enough getopt to handle long options (see the man page)
getopt -T &>/dev/null
if [[ $? -ne 4 ]]; then echo "Getopt is too old!" >&2 ; exit 1 ; fi

# read the options
TEMP=`getopt -o '' --long 'sample:,read1:,read2:,genome:,work:,outdir:' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        --sample) sample_name=$2 ; shift 2;;
        --read1) read1=$2        ; shift 2;;
        --read2) read2=$2        ; shift 2;;
        --genome) genome=$2      ; shift 2;;
        --work) work=$2          ; shift 2;;
        --outdir) outdir=$2      ; shift 2;;
        --) shift ; break ;;
        *)
        echo ""
        echo "Error in given Parameters. Undefined: "
        echo $*
        echo ""
        exit 1
    esac
done

# Print information about the input parameters.
echo "Mapping reads for sample $sample_name"
echo "Path to file with mate 1 reads: $read1"
echo "Path to file with mate 2 reads: $read2"

echo "Mapping reads to: $genome"
echo "Setting working directory to $work"
echo "Using trim_galore arguments: $trim_galore_args"
echo "Using Bismark arguments: $bismark_args."

# Trimming with trim galore!
trim_galore_dir=$work/trimgalore
echo ""
echo "Trimming reads, and saving output to: $trim_galore_dir"
echo "Running fastqc on the trimmed reads, and saving to $work/fastqc"
mkdir -p trim_galore_dir
mkdir -p $work/fastqc

trim_galore \
    --gzip \
    --paired \
    -o $trim_galore_dir \
    --basename $sample_name \
    --fastqc \
    --fastqc_args "--outdir $work/fastqc" \
    $trim_galore_args \
    $read1 $read2

echo "Finished trimming reads."
date
echo ""


# Align reads to the genome
bismark_dir=$work/bismark/$sample_name
mkdir -p $bismark_dir
echo "Aligning reads using Bismark and saving the output to"
echo "$work/bismark"

bismark $(dirname $genome) \
    -1 $trim_galore_dir/$sample_name'_val_1.fq.gz' \
    -2 $trim_galore_dir/$sample_name'_val_2.fq.gz' \
    ${bismark_args} \
    -o $bismark_dir

echo "Finished aligning reads."
date


# sort the newly aligned BAM file.
sorted_by_name=$bismark_dir/$sample_name.sortedByName.bam
mkdir -p $(dirname $sorted_by_name)
echo "Sorting the newly aligned BAM file by read name, and saving as:"
echo $sorted_by_name

samtools sort \
    -n \
    -@ $(nproc) \
    -o $sorted_by_name \
    $bismark_dir/$sample_name'_val_1_bismark_bt2_pe.bam' # Need to fill this in!
echo "Finished sorting."
date


# Deduplicate BAM file
mkdir -p $(dirname $sorted_by_name)
echo Deduplicating the sorted BAM file and saving as $sample_name.sortedByName.deduplicated.bam

deduplicate_bismark \
    --bam \
    --paired \
    --outfile $sample_name \
    -output_dir $bismark_dir \
    $sorted_by_name

echo "Finished deduplication."
date


# Sort and index the deduplicated file by position
sorted_by_pos=$bismark_dir/$sample_name.sortedByPos.bam
mkdir -p $(dirname $sorted_by_pos)
echo "Sorting the deduplicated BAM file by position, and saving as:"
echo $sorted_by_pos
samtools sort \
    -@ $(nproc) \
    -o $sorted_by_pos \
    $bismark_dir/$sample_name.deduplicated.bam
echo "Finished sorting. Indexing..."
samtools index  $sorted_by_pos
echo "Finished indexing."
date


# extract methylation calls
echo "Running bismark_methylation_extractor, and saving to:"
echo $bismark_dir/reports
genome_folder=$(dirname $genome)
meth_extractor=$scratch
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


# Stage out the files to keep
mkdir -p $outdir/logs
cp $scratch/$sample_name.deduplication_report.txt $outdir/logs
cp $scratch/$sample_name_val_1_bismark_bt2_PE_report.txt $output_dir/logs

mkdir -p $outdir/reports
mv $scratch/reports/$sample_name.deduplicated.CX_report.txt.gz $outdir/reports/$sample_name.CX_report.txt.gz

mkdir -p $outdir/sorted
cp $scratch/$sample_name.sortedByPos.bam $outdir/sorted
cp $scratch/$sample_name.sortedByPos.bai $outdir/sorted