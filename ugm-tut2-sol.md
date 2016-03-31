---
layout: default
permalink: ugm/tut2/solution/
---

[Go back](/ugm/tut2/)

# Tutorial 2 - Solution

The output file for translated clonotype frequencies is specified using `-oa`. Hence, the first configuration corresponds to the following call of **IMSEQ**:

~~~Plaintext
$ imseq -ref Homo.Sapiens.TRB.fa -j 4 -oa DEFAULT examples/data/example_quality_bias.fq.gz
~~~

The additionally specified option `-j 4` is optional and causes **IMSEQ** to use four parallel treads. You can adjust this depending on the hardware available.

To reduce the quality threshold, the option `-mq` has to be specified:

~~~Plaintext
$ imseq -ref Homo.Sapiens.TRB.fa -j 4 -mq 10 -oa MQ10 examples/data/example_quality_bias.fq.gz
~~~

and to additionally enable simple and quality clustering, use `-sc` and `-qc`:

~~~Plaintext
$ imseq -ref Homo.Sapiens.TRB.fa -j 4 -mq 10 -sc -qc -oa MQ10SCQC examples/data/example_quality_bias.fq.gz
~~~

To plot the top 10 clonotypes of `MQ10SCQC` and compare their frequencies to those of the other generated files, you can use the provided R script:

~~~Plaintext
./examples/top10PlotPDF.R MQ10SCQC MQ10 DEFAULT plot.pdf
~~~

This should result in the following plot:

<div style="text-align:center; margin-bottom:20px">
<img alt="Tutorial 2 Result" src="/images/tut2plot.png"/>
</div>

The following clonotype counts can be observed:

~~~Plaintext
$ wc -l DEFAULT MQ10 MQ10SCQC
     283 DEFAULT
     437 MQ10
     241 MQ10SCQC
~~~

With a strict quality threshold of 30, 283 clonotypes are detected. Relaxing that threshold allows for more low quality reads to be analyzed, resulting in a more diverse repertoire of 437 clonotypes. Correcting for errors in that repertoire reduces the size to 241 clonotypes.
