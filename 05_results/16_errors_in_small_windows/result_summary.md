# Can we distinguish CG-only, TE-like and non-methylation?

**Date:** 9th August 2023
**Author:** Tom Ellis

## Background

Result 15 found that conversion error rates across 1000bp windows of the 
chloroplast vary quite a bit.
There also seems to be periodicity in the error rates:
    - In result 11 I found autocorrelation over ~20bp
    - Greg found cycles of around 150bp. This is very similar to the scale at 
        which tn5 cuts
    - The 1kb windows might show cycles as at larger scales as well.

We suspect that errors are caused by tn5 making single-strand nicks roughly 
every 150bp.
That means 1kb is too large a scale to look at.
Here I want to look at 150bp windows and reestimate the distribution of errors

## What did you do?

- `01_get_methylation_counts.py` is a generic script calling `methylation_in_windows`
- `02_submit_methylation_counts.sh` runs the python script on four Col0 samples.
- `03_multiple_chloroplasts.R` plots and estimates beta parameters for the four
  error distributions.
- `04_inferring_TE_status.qmd` repeats a previous analysis from result 13 to
  call methylation status (unmethylated, CG-only, TE-like) using updated beta
  parameters

## Main conclusion
Beta distributions are wider than using 1kb windows, but distributions look similar.
Mean across four samples: a=8.529238, b=111.4319

TE calls on different samples look very similar than using 1kb windows.
Gene calls look a bit more like what I expected from Zhang et al
One window on the chloroplast looks TE-methylated.

## Caveats

In results 13 a large part of the chloroplast had no reads mapping.
Here the whole chloroplast is covered. I don't know why.

## Follow-up
