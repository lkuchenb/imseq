---
layout: default
permalink: ugm/tut3/solution/
---

[Go back](/ugm/tut3/)

# Tutorial 3 - Solution

To perform the regular analysis without barcode correction, the reads can either be preprocessed, i.e. using `fastx_trimmer` from the [FASTX Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) or by specifing `-bcl 10` (causing **IMSEQ** to interpret the first 10 bases as barcode) in combination with `-ber 0` (not allowing for any correction):

~~~Plaintext
$ imseq -ref Homo.Sapiens.TRB.fa -j 4 -bcl 10 -ber 0 -oa BC-default examples/data/example_barcode_correction.fq.gz
~~~

The optional flag `-j 4` enables **IMSEQ** to perform the barcode based correction and the clonotyping using 4 threads.

Without `-ber 0` the default maximum error rate of `0.05` is used. To additionally write the barcode count based clonotype counts, use `-oab`:

~~~Plaintext
$ imseq -ref Homo.Sapiens.TRB.fa -j 4 -bcl 10 -oa BC-corr -oab BC-corr-dist examples/data/example_barcode_correction.fq.gz
~~~

You can see that the barcode based correction reduced the number of clonotypes from 615 to 575:

~~~Plaintext
$ wc -l BC-*
     575 BC-corr
     575 BC-corr-dist
     615 BC-default
~~~

Using the provided R script, you can visualize the distribution of the top 10 clonotypes:

~~~Plaintext
$ ./examples/top10PlotPDF.R BC-default BC-corr BC-corr-dist plot.pdf
~~~

This should result in the following plot:

<div style="text-align:center; margin-bottom:20px">
<img alt="Tutorial 3 Result" src="/images/tut3plot.png"/>
</div>

One can make two observations:

 * Using the PCR barcodes to correct for PCR and sequencing errors slightly boosts the head of the distribution, since the low frequency erroneous clonotypes in the tail of the distribution get corrected.
 * Deriving the distribution from the unique barcode counts per clonotype flattens it, indicating that the highly dominant clonotypes partially show up that high due to preferential PCR amplification. Even in the top 10, some rank shifts occur compared to the uncorrected distribution.
