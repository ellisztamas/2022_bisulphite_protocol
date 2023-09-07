# Can we distinguish CG-only, TE-like and non-methylation?

**Date:** 24th July 2023
**Author:** Tom Ellis

## Background

Bob schmitz makes a solid case that DNA exists in a continuum between
unmethylated DNA, CG-only methylation and TE-like methylation (mC in all three
contexts).
I described a framework to work out which is which in 06_reports and ran
simulations.
Here I want to confirm that it works on real data by comparing cases where we
probably know the right answer by looking at:
    - chloroplast (unmethylated)
    - gene bodies (these should be a mixture of stuff tbh...)
    - TEs known to be methylated by the RdDM and CMT2 pathways
## What did you do?

- I looked at a single genome of Col-0 sequenced at 30x coverage
- `01_get_methylation_counts.py` gets methylation levels on the TEs, genes, and
  in 1000-bp windows across the genome
- `02_submit_methylation` runs that as a job array.
- `03_inferring_TE_status.R` calls methylation status calculates the probability
  that observed methylated read counts are drawn from the (beta) distribution
  of errors, or from more than that (the CDF of the same beta distribution)

## Main conclusion
When we allow for residual errors we can recover unmethylated status across a 
chloroplast reliably.

Likewise, most RdDM and CMT2 TEs come back as TE-methylated.
The ones that don't have (un)methylated read counts that make sense, so this is
either real or something to do with mapping.

46% of genes appear unmethylated; 38% CG methylated15% TE-like
Zhang et al (2020) suggest it should be more like 80%, 15, 5.
So, probably some false calls, although their data may not to be perfect either.

## Caveats

Likely that methylation still overestimated, and calls are not as conservative,
as other papers

## Follow-up
Comparing the PDF and CDF of the beta distribution basically asks if data are
consistent with (1) the error distribution or (2) more than that.
I can probably do better by using more informative priors about the expected 
methylation levels for each context.

Check smaller window sizes.

Do an HMM!