# Title Separating original and synthesised strands

**Date:** 21st December 2022
**Author:** Tom Ellis

## Background
In BS libraries one distinguishes between the original top and bottom (OT, OB) 
strands, the complementary strands synthesised during PCR (CTOT, CTOB).
I wanted to see if conversion errors were different on the original vs 
complementary strands.

## What did you do?
The bioinformatician on the protocol paper wrote a script to separate these
reads. I ran this, then ran the nextflow pipeline on the original and synethsised
reads separately

## Main conclusion
You get mistakes on the oirginal and complementary strands.
Errors do not result from funkiness during PCR.

## Caveats

## Follow-up
