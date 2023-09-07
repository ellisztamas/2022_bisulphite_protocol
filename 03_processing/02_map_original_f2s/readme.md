Map reads for the original samples corresponding to those in the mix plate.

These are the original samples from Rahul Pisupati's temperature experiment 
sequenced using the 'old' protocl using unmethylated dNTPs that were subsequetly
redone with the 'new' protocol.

I have not repeated the creation of the fastq files from raw bams or the sample
sheet because this was done in another project already. See this directory for
details of how that was done:
/groups/nordborg/projects/epiclines/007.hmm/03_processing/
Instead, it just subsets the sample sheet from that directory and remaps the
reads using the methylseq pipeline to the TAIR10 reference genome (cf the
`007.hmm` folder, which maps reads to the parental genomes of the cross.)