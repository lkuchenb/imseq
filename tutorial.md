---
layout: default
date:   2015-02-04 15:53:59
permalink: tutorial/
---

# Tutorial

**IMSEQ** comes with two example input files and the reference segment files for the human T cell receptor alpha (TRA) and beta (TRB) chain genes. 

 * TOC
{:toc}

## Example 1: Simulated error-free data

The first example is an error-free simulated dataset of 50 full gene sequences. The data can be found in ```examples/example_input_sim.fa```. To process the data, **IMSEQ** can be invoked as follows:

    $ imseq -r -ref ../Homo.Sapiens.TRB.fa -o output.tsv example_input_sim.fa

The switch ```-r``` has to be specified since this is a simulated sequence and not a real reverse complement read from a targeted enrichment PCR + sequencing experiment. The ```-ref``` option is used to tell **IMSEQ** where to find the FASTA file with the [correctly formatted](../fastaFormat) V and J segment sequences. With ```-o``` we instruct **IMSEQ** to write an output file with detailed information about every successfully processed read. The last argument passed to **IMSEQ** is the input file. 

The file ```output.tsv``` now contains TAB separated information about every input read:

    cdrBegin  cdrEnd  leftMatches  leftErrPos  leftMatchLen  rightMatches  rightErrPos  rightMatchLen  cdrNucSeq                                cdrAASeq
    270       303     V6-2V6-3     -           273           J2-4          -            31             TGTGCCAGCATTTCCAAAGACATTCAGTACTTC        CASISKDIQYF
    270       297     V3-1         -           273           J2-1          -            31             TGTGCCTGCTACAATGAGCAGTTCTTC              CACYNEQFF
    273       303     V7-4         -           276           J2-1          -            31             TGTGCCGGCTGCTTAGATGAGCAGTTCTTC           CAGCLDEQFF
    273       312     V23-1        -           276           J2-6          -            31             TGCGCCAGCAGGCAATCTGAGGCCAACGTCCTGACTTTC  CASRQSEANVLTF
    257       284     V11-3        -           260           J1-4          -            31             TGTGCCACCAGTGAAAAACTGTTTTTT              CATSEKLFF

<p style="margin-top:50px"></p>

## Example 2: Low quality real world data

The second example file contained in the **IMSIM** downloads is ```examples/quality_bias_sample.fq.gz```. It originates from a low quality sequencing run and illustrates the impact of low quality reads on the clonotype distribution when they are filtered out. We will process the data twice, once with a quality threshold of 30 and no posterior clustering and once with a quality threshold of 10 with posterior clustering:

    $ imseq -ref ../Homo.Sapiens.TRB.fa -j 8 -mq 10 -qc -oa MQ10QC quality_bias_sample.fq.gz
    $ imseq -ref ../Homo.Sapiens.TRB.fa -j 8 -mq 30 -oa MQ30 quality_bias_sample.fq.gz

This time, we do not write a full per-read output file but aminoacid based clonotype counts (```-oa```). **IMSEQ** shall also use 8 parallel threads (```-j 8```). A clonotype count output file look like this:

    $ sort -rnk2 MQ30 | head
    V20-1:CSARDELAANYEQYF:J2-7            630
    V6-5:CASSYSTGVIDGYTF:J1-2             592
    V6-5:CASSFQTGTGVYGYTF:J1-2            286
    V12-3V12-4:CASSTQGYEQYF:J2-7          131
    V10-3:CAVRLAGETDTQYF:J2-3             119
    V27:CASSPTSGAPGEQFF:J2-1J2-5          90
    V2:CASMAGTYNEQFF:J2-1J2-5             87
    V12-3V12-4:CASSAQETQYF:J2-5           87
    V6-1:CASRGYNEQFF:J2-1J2-5             80
    V12-3V12-4:CASSAQETQFF:J2-1           72

It contains one row for each clonotype associated with its count, separated by a TAB character. This data can be used to examine the two different clonotype distributions generated with the two different parameter sets. We provide a small R script in the examples directory, ```examples/top10PlotPDF.R```, which visualizes the frequencies of the top 10 clonotypes in a bar diagram. It is invoked as follows:

    $ ./top10PlotPDF.R MQ10QC MQ30 MQ.pdf

The output file, ```MQ.pdf```, should contain the following plot:

![Tutorial Example 2 Result](../images/MQ30vsMQ10example.png)

The example clearly shows that a simple quality filtering approach can drastically change the clonotype distribution and indicates that erroneous reads should rather be rescued and corrected than discarded.
