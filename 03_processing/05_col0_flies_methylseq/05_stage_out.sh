# Stage the cytosine coverage reports, Bismark alignment reports and multiqc
# report out of scratch-cbe.

scratch=/scratch-cbe/users/$(whoami)/05_col0_flies/methylseq
target_dir=03_processing/04_col0_flies_methylseq/output
mkdir -p $target_dir

stage $scratch/bismark/sorted $target_dir
stage $scratch/bismark/coverage2cytosine/reports $target_dir
stage $scratch/bismark/alignments/logs $target_dir
cp $scratch/multiqc/bismark/multiqc_report.html $target_dir
mv $target_dir/multiqc_report.html $target_dir