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

## Main conclusion

## Caveats

## Follow-up
