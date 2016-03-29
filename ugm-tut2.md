---
layout: default
permalink: ugm/tut2/
---

# Tutorial 2

IMSEQ can collapse clonotype clusters based on quality scores or simply based on user defined distance thresholds. To illustrate the effect of this correction, we will now analyse a real world sample once with and once without the correction enabled. The sample data can be found in the file `examples/data/example_quality_bias.fq.gz`, which contains 10,000 reads suffering from unusually low quality.

The relevant options can be found in the corresponding sections of the help output (except):

~~~Plaintext
Quality control:
  -mq, --min-qual NUM
        Minimum average read phred score. In paired end mode, this is
        applied to both reads. See '-sfb'. In range [0..60]. Default: 30.

Postprocessing / Clustering:
  -qc, --qual-clust
        Enable quality score based clustering.
  -sc, --simple-clust
        Enable simple distance-based clustering
  -qcme, --max-err-hq NUM
        Maximum edit-distance for two clusters to be clustered without low
        quality correlation In range [0..inf]. Default: 4.
  -qcsd, --min-sd-diff NUM
        How many standard deviations must an error positions quality value
        be below the mean to be considered for clustering. Default: 1.
  -scme, --max-err-lq NUM
        Maximum edit-distance for two clusters to be potentially clustered
        without low quality correlation In range [0..inf]. Default: 2.
  -mcr, --max-clust-ratio NUM
        Maximum abundance ratio for two clonotypes to be clustered In range
        [0..1]. Default: 1.
~~~

As stated in the help output, by default IMSEQ discards reads with an average quality score below 30.

## Assignment

1. Analyse the supplied dataset with **IMSEQ** in the following three configurations. Write a file with translated clonotype frequencies for each analysis.
    1. Use the default parameters
    1. Reduce the average quality threshold to 10
    1. Reduce the average quality threshold to 10 and enable both simple and quality clustering
1. Plot the distribution of the top 10 clonotypes using the provided R script `top10PlotPDF.R`. Supply the lastly generated file first.
1. How many clonotypes were detected with each configuration?

Click [here](/ugm/tut2/solution/) to see the solution.
