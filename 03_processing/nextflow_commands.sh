# Load nextflow modules
module load build-env/f2022
module load nextflow/22.10.0

# Tell trim galore explicitly to trim the first 15 bp from each read
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/00_temp/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/wgbs \
--trim_galore_args="--clip_r1 15 --clip_r2 15" \
-with-conda = true

# Tell trim galore explicitly to trim the first 15 bp from each read
# Remove flag `--non-directional` in call to bismark
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/00_temp/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/directional \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --three_prime_clip_R1 15 --three_prime_clip_R2 15" \
--bismark_args="--score_min L,0,0.5" \
-with-conda = true

# Tell trim galore explicitly to trim the first 15 bp from each read
# and also that reads are paired
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/00_temp/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/wgbs_paired \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --paired" \
-with-conda = true

# Tell trim galore explicitly to trim the first 15 bp from each read
# and to keep only reads >100bp
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/00_temp/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/wgbs_min_100bp \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --paired --length 100" \
-with-conda = true

# Tell trim galore explicitly to trim the first 15 bp from each read 
# in both the 5- and 3-prime directions, keeping reads > 100bp
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/00_temp/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/trim_r1_min_100bp \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --paired --length 100 --three_prime_clip_R1 15 --three_prime_clip_R2 15" \
-with-conda = true

# Tell trim galore explicitly to trim the first 15 bp from each read 
# and Bismark to treat dovetailing reads as non-concordant
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/00_temp/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/no_dovetail \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --paired" \
--bismark_args="--no_dovetail" \
-with-conda = true


# Tell trim galore explicitly to trim the first 15 bp from each read 
# in both the 5- and 3-prime directions
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/00_temp/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/three_prime_15 \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --three_prime_clip_R1 15 --three_prime_clip_R2 15" \
-with-conda = true 

# Tell trim galore explicitly to trim the first 15 bp from each read 
# in both the 5- and 3-prime directions
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/00_temp/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/trim_50 \
--trim_galore_args="--clip_r1 50 --clip_r2 50 --three_prime_clip_R1 50 --three_prime_clip_R2 50" \
-with-conda = true 

# Tell trim galore explicitly to trim the first 15 bp from each read 
# in both the 5-prime directions, and 15-bp from read 1 in the 3-prime direction
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/00_temp/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/trim_r1 \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --three_prime_clip_R1 15" \
-with-conda = true 

# Tell trim galore explicitly to trim the first 15 bp from each read 
# in both the 5-prime directions, and 15-bp from read 2 in the 3-prime direction
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/00_temp/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/trim_r2 \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --three_prime_clip_R2 15" \
-with-conda = true 


# Run this on pentuple mutants that shouldn't have any DNA at all
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "/scratch-cbe/users/thomas.ellis/pentuple_mutant/000000000-GDJK9_0_R14404_20221025/demultiplexed/211006/*_R{1,2}_*.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/mutants \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --three_prime_clip_R1 15 --three_prime_clip_R2 15" \
--bismark_args="--score_min L,0,0.5 --non_directional" \
-with-conda = true 


# An example from Rahul's data
# Doesn't work.
nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "01_data/01_old_protocol/HHHY7DSX2_2#169198_ACTCGCTAAAGGAGTA.fastq" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 04_output/old_protocol \
--trim_galore_args="--clip_r1 15 --three_prime_clip_R1 15" \
--bismark_args="--score_min L,0,0.5" \
--single_end = true \
-with-conda = true
# -w /scratch-cbe/users/thomas.ellis