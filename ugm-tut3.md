---
layout: default
permalink: ugm/tut3/
---

# Tutorial 3

The file `example_barcode_correction.fq.gz` provided with **IMSEQ** contains 20,000 reads from an experiment, where 10 random nucleotides were incorporated into the gene copies generated when enriching the target genes. The provided sequences are VDJ reads starting with the 10 random bases.

To instruct **IMSEQ** to perform barcode based error corrections, the `-bcl` flag has to be used to set a barcode length other than 0 (default). The PCR barcode related flags are shown below:

~~~Plaintext
Barcoding:
  -bvdj, --barcode-vdj
        In paired end mode: Read the barcode from the VDJ read instead of
        the V read.
  -bse, --bcseq-max-err NUM
        Maximum number of errors allowed in the barcode sequence In range
        [0..inf]. Default: 1.
  -bmq, --bc-min-qual NUM
        Minimum per base quality in molecular barcode region In range
        [0..60]. Default: 30.
  -bcl, --barcode-length NUM
        Length of random barcode at the beginning of the read. A value of
        '0' disables barcode based correction. In range [0..inf]. Default:
        0.
  -ber, --barcode-err-rate NUM
        Maximum error rate between reads in order to be merged based on
        barcode sequence In range [0..1]. Default: 0.05.
  -bfr, --barcode-freq-rate NUM
        Inclusive maximum frequency ratio between smaller and larger cluster
        during barcode clustering In range [0..1]. Default: 0.2.
  -bst, --barcode-stats FILE
        Path to barcode stats output file. If empty, no file is written.
        Ignored if -bcl is 0.
  -oab, --out-amino-bc FILE
        Output file path for translated clonotypes with barcode corrected
        counts. Ignored if -bcl is 0.
  -onb, --out-nuc-bc FILE
        Output file path for untranslated clonotypes with barcode corrected
        counts. Ignored if -bcl is 0.
~~~

When enabled, the PCR barcode error correction corrects PCR and sequencing errors in the observed input sequences according to the specified parameters `-ber` and `-bfr`.

Additionally, an extra output file format is available through the flags `-oab` and `-onb`. It corresponds to `-oa` and `-on`, i.e. the (translated or nucleotide based) clonotype counts, but does not contain the number of reads corresponding to each clonotype but the number of unique barcodes. These can be helpful in order to avoid amplification biases altering the observed distributions, i.e. if certain V or J segments specific primers perform better than others.

## Assignments

1. Analyse the provided example dataset `example_barcode_correction.fq.gz` in two configurations:
    1. Using the default analysis mode, i.e. without barcode correction. Write the translated clonotype counts to a file. *Note that the barcodes disrupt the analysis if IMSEQ interprets them as a valid part of the gene - you can either preprocess the file by removing them with a tool of your choice or instruct **IMSEQ** to read the barcode but correct for no errors.*
    1. Using the barcode analysis mode. In this case, write both the uncorrected and the corrected translated clonotype counts to a file.
1. How many unique clonotypes are detected with and without barcode based error correction?
1. Compare the distributions of the top 10 clonotypes of the three produced files using the provided script `top10PlotPDF.R`.

Click [here](/ugm/tut3/solution/) to see the solution.
