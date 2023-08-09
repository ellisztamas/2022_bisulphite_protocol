# New bisulphite-sequencing protocol, 2022

In 2022 Viktoria spent a lot of time formulating and optimising a new
bisulphite-sequencing protocol to save time and money. Concerns were raised 
that conversion rates came back as being suspiciously low (~5%).

Note: I (Tom Ellis) am writing this in November 2022, having been on paternity
leave for most of the year while the experimental procedures were done, so
I am coming to this late, and being not terribly knoweldgeable about lab and 
bioinformatics details in any case. As such this text might contain errors.

## Why could conversion rates be low?

A summary of the bisulphite protocol:

1. DNA is sheared into fragments, ideally of 200bp.
2. Adapters are ligated onto the 5\`- and 3\` ends. This entails a gap-
    repair step, incorporating dNTPS. The new protocol uses methylated dCTPs
3. Bisulphite conversion. The dCTPs aren't converted, saving a step.
4. Sequencing with PE100, which in the ideal cases sequences 100bp from the 5- 
    and 3-prime ends, to cover the whole 200bp fragment
5. Align reads to the genome, including the chloroplast genome. Since this 
    should be unmethylated, methylation on the chloroplasts give an estimate of
    the genome-wide conversion rate.
6. ????
7. Publish.

Where can errors happen:

- When fragments are much shorter than 200bp the forward sequence can include
    the 3\` adapter. 
    - Signal: GC content increases continuously towards the 3\`-end of the reads
    - Solution: ensure the pipeline looks for and removes adapter sequences.
- The gap repair step incorporates random nucleotides
    - Signal: GC content at the 5\`-end of reads is unstable (lines go up and 
    down).
    - Solution: remove the first 14 base pairs from each read.
- DNA fragments break somewhere in the middle, because BS treatment is brutal, 
    and are repaired, incorporating random dNTPs.
    - Signal: GC content is uneven across the read
    - Solution: cry, then do it again.

Data from test runs of the new protocol done throughout 2022 came back with 
lower apparent conversion rates than recent data using the previous protocol
(specifically, Rahul's data on F2 crosses grown at 4°C and 16°C).

## How to navigate this folder

### Data

See the readme file in `01_data` for details of different datasets.

### Mapping reads

The primary effort was to fix the non-conversion issue by playing with the 
bioinformatics settings.

For the commands used see `03_processing/nextflow_commands.sh`.
Output is in `04_output`.

I used a fork of the `nextflow_pipelines` repo from Rahul, who ultimately forked
it from Felix Krüger and Simon Andrews. 
See https://github.com/ellisztamas/nextflow_pipelines
I used the `nf_bisulphite_WGBS` module, which runs:

1. Trim Galore! to trim reads
2. Fastqc on reads before and after trimming
3. Bismark to align and call methylation
4. Samtools to sort stuff. I didn't totally follow what it's doing
5. Multiqc to report the output

So far I used a custom hacky script to get conversion rates manually.
See `03_processing/conversion_rates_manually.py`.

### Other results

For more involved results that went beyond just fiddling with the trimming
settings see the folders in `05_results`. Each folder has a `result_summary.md`
file that summarises what I did and the conclusions.

See also 06_reports for more involved thinking about how to deal with the
low conversion rates.

## Summary of conclusions

- For 5 samples from the spike-in done on plate 2021-015 I tested various ways to
    trim reads from the 5 and 3 prime ends of reads:
    - Before trimming, you can see the effect of random nucleotides from the gap
        repair step, and adpater content.
    - Trimming 5 prime ends removes this, but gives conversion rates ~4%.
    - Trimming 3 prime ends reduces conversion to ~2%, which is still not enough.
    - Conversion rates are similar for the chloroplast and Lambda vector, suggesting 
        that DNA complexity is not the issue.
- Pentuple mutant:
    - All samples showing conversion rates of 10e-2
    - This includes the Col-0 control
    - This did not depend on transposase concentration.
- There is an odd periodicity to the mistakes.
    - I found autocorrelation in the position of methylated cytosines on the
        chloroplast over 20 bp
    - Greg found clear periodicity over ~150 bp.
    - This is also true in the lambda vector, so it is unlikely to be residual 
        chromatin packing.
- We asked Bob Schmitz his opinion:
    - in addition to making double strand breaks, the tn5 transposase could be
        making single-stand nicks all over the place, that are then repaired
        with methylated dNTPs.
    - The 150bp business might also have something to do with the size of tn5,
        and that it jumps along the chromosome at 150bp intervals.
