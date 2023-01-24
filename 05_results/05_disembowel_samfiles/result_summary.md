# Title Disembowel SAM files

**Date:** 18th January 2023
**Author:** Tom Ellis

## Background

I suspect that high conversion rates arise because with a low probability strands
get damaged and are repaired in the presence of methylated dNTPs. That causes 
the read to have a whole string of methylated cytosines that are *next to one
another*.

I want to get a distribution of numbers of methylated cytosines on each read and
quantify how often those methylated cytosines are found as a cluster.

## What did you do?

As a toy example, I took the reads from the chloroplast of one sample 
(generated with `extract_chloroplast.sh`) and went through how to take the file
apart in detail, based on tips from Yoav Voichek. See `disembowel_samfiles.ipynb`)

Then I systematically went through 4 aligned BAM files and counted how many
methylated cytosines were on each, and how often those cytosines cluster into a
contiguous string.
- `meth_counts.py` is a Python script to parse a BAM, calling library module
    `BismarkSam.py`.
- `meth_counts.sh` runs the script on each BAM file and exports to `meth_counts`
- `meth_counts.Rmd` plots the results.

## Main conclusion
60-70% of reads on the chloroplast have cytosines clustered together, which is 
not as many as I expected.
They tend to be the very longest reads.

I also checked IGV again visually, and it does seem like getting rid of
clustered reads should solve most of the problem.

## Caveats
I tried to do a comparison with Rahul's folder, but it seems the anatomy of
those SAM files is different.

## Follow-up
Compare with the old protocol.
Go through and pull out read *pairs*, because I think one bad read will pull
out a problematic pair that would otherwise be easy to detect.