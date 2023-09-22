# Stage the cytosine coverage reports, Bismark alignment reports and multiqc
# report out of scratch-cbe.

scratch=/scratch-cbe/users/$(whoami)/original_f2s/methylseq/
target_dir=03_processing/03_map_original_f2s/output
mkdir -p $target_dir

stage $scratch/bismark/sorted $target_dir
stage $scratch/bismark/coverage2cytosine/reports $target_dir
stage $scratch/bismark/alignments/logs $target_dir
cp $scratch/multiqc/bismark/multiqc_report.html $target_dir
mv $target_dir/multiqc_report.html $target_dir/
