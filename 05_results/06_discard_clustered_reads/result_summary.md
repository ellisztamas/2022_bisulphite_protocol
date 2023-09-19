# Discard reads with clustered methylated cytosines 

**Date:** 24th January 2023
**Author:** Tom Ellis

## Background
We suspect that methylated cytosines are being incorporated into reads when they
are nicked by Bst. I previously counted how often reads have cluster of 
contiguous cytosines.

Here I want to go through each pair of reads in a BAM file, check if one of the 
pair has a string of contiguous mC and keep if neither does.

## What did you do?
- `discard_reads.sh` goes through every BAM file in `04_output/three_prime_15/bismark_dedup`
    and:
    - converts to SAM
    - runs `discard_clustered_reads.py` to go through (pairs of) reads in each 
    sam file and check whether either mate has a strong of contiguous methylated
    cytosines.
    - Runs the trimmed sam files through `bismark_methylation_extractor`.
    - runs `02_library/calculate_conversion.py` on resulting allc files
- See `conversion_rates.csv` for the output
- Also applied the conversion-rate script to the output `04_output/three_prime_15/` and 
    `04_output/wgbs_paired/`.
- I visually compared sorted bam files for cleaned vs discarded reads, and they
    show the pattern you would expect.

## Main conclusion

Trimming reads with contiguous cytosines brings non-conversion rates down from
0.019-0.035 to 0.006-0.011.

## Caveats
There is quite a bit of variation between samples, and I only ran it on five files.

Methylation levels on the other chromosomes goes *up* from ~8% to ~25%

Need to double check that the check for whether mate pairs match makes sense!
I used the read tag, but check the bitwise code

## Follow-up
Compare with data from the old protocol