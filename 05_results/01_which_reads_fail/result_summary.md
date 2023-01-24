# Result name

**Date:** 19th December 2022
**Author:** Tom Ellis

## Background
We want to know whether unusual conversion rates are due to a small number of reads doing something funky

## What did you do?

For the allc files in `04_output/wgbs_paired`, pull out the cytosines on the
Lambda_NEB and chloroplasts only and filter for positions with at least one 
methylated cytosine.

## Main conclusion

* On the Lambda vector it is clear that most methylated positions are clustered
to within the size of a read. Because coverage is generally <1 per cytosine,
they must be from the same read.
* Coverage is much higher on the chloroplast so the pattern is not as crystal
clear, but still present.

## Caveats

## Follow-up

Find a way to dig deeper into which reads are messed up.