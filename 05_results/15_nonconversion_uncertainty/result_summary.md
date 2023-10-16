# Can we distinguish CG-only, TE-like and non-methylation?

**Date:** 24th July 2023
**Author:** Tom Ellis

## Background

It seems that if we have a good estimate of the conversion error rates we can
accurately correct methylation estimates.
We can estimate conversion rates from reads mapping to the chloroplasts.

However, is that a good estimate for the rest of the genome?
Does the error rate vary across the genome?
Can we fix that?

## What did you do?

- `01_get_methylation_counts.py` gets methylation levels on the chloroplasts,
    CMT2- and RdDM-targetted TEs and genes for 4 cytosine report files from
    Col-0 at high coverage.
- `02_submit_methylation` runs that as a job array.
- `03_chloroplast.R` pokes about in them and checks whether error rates are 
    uniform, and whether we can correct for that
- `02_library/R/binomial_with_errors.R` gives functions to correct for errors

## Main conclusion

- There is substantial variation in error rate across a chloroplast above what 
    would be expected from binomial sampling
- Patterns are repeatable between technical replicates
- Patterns are not explained by annotated features

## Caveats

Chloroplasts don't have any TEs

## Follow-up

Repeat for a different biological replicate
Check windows at smaller scales.