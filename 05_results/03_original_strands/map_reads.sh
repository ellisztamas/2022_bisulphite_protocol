module load build-env/f2022
module load nextflow/22.10.0

nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "05_results/03_original_strands/read{1,2}.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 05_results/03_original_strands/pipeline_output \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --three_prime_clip_R1 15 --three_prime_clip_R2 15" \
-with-conda = true

nextflow run 02_library/nextflow_pipelines/nf_bisulfite_WGBS \
--input_files "05_results/03_original_strands/read{1,2}_unassigned.fastq.gz" \
--fasta 01_data/03_reference_genome/TAIR10_wholeGenome_withVectors.fasta \
--outdir 05_results/03_original_strands/pipeline_output_unassigned \
--trim_galore_args="--clip_r1 15 --clip_r2 15 --three_prime_clip_R1 15 --three_prime_clip_R2 15" \
-with-conda = true

