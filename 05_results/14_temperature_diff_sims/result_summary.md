# Can we distinguish mean methylation at two temperatures?

**Date:** 27th July 2023
**Author:** Tom Ellis

## Background

Laura wants to characterise methylation levels of TEs in plants grown at cold and warm temperatures.
This follows on from Manu's (2016) paper. Looking at that paper:

* CHH is the only sequence where we expect to see detectable changes.
* figure 1 suggests genome-wide CHH seems to go from about 18% at 10°C to 22% at 16°C
* figure 4 suggests that CHH on (variable) TEs varies between about 12% and 22% across the TE body, but is a consistently few percent higher at the warmer temperature.

It is immediately obvious that how well we can tell apart methylation depends on the difference between the temperatures, the absolute magnitude of methylation (things closer to zero or one will be harder), and number of trials (coverage times number of cytosines).

## What did you do?

- `02_library/R/binomial_with_errors.R` gives R functions to simulate a binomial
    process with errors:
        - generate random draws
        - estimate likelihood of data
        - return maximum likelihood value
- `temp_diff_sims.R` simulates methylation
    at temperature a, and at temperature b with is *d* higher than a. This is
    done for different values of:
        - d
        - mean at temperature a
        - error rates
        - number of trials (number of cytosines * coverage)
    This uses only a point estimate for error rates
- `plot_diff_sims.R` plot results of the above.
- `average_over_TEs.R` is a simulation to ask how many TEs one would have to 
    aggregate to reliably estimate a mean methylation

## Main conclusion

Binomial sampling variance is a larger issue that errors.
You need to average over >100 TEs to reliably estimate a mean

## Caveats

Assumes the TEs you aggregate have the same true methylation level.

## Follow-up