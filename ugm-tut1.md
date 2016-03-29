---
layout: default
permalink: ugm/tut1/
---

# Tutorial 1: Basic usage of IMSEQ

The basic usage of **IMSEQ** is displayed on the command line when **IMSEQ** is invoked without any arguments:

~~~Plaintext
imseq
=====
    imseq -ref <segment reference> [OPTIONS] <VDJ reads>
    imseq -ref <segment reference> [OPTIONS] <V reads> <VDJ reads>
    Try 'imseq --help' for more information.
~~~

The `-ref` switch has to be specified and point to a file with the V and J gene segments of the corresponding gene, [formatted in a special way](/fastaFormat/). **IMSEQ** comes with reference files for the human genes TRA, TRB, IGL, IGH and IGK. Additionally, at least one of the output switches `-o`, `-oa` and <nobr markdown="span">`-on`</nobr> has to be specified.

The file `examples/data/example_sim.fa` contains 50 simulated full length human TRB gene sequences. The corresponding reference segment file is `Homo.Sapiens.TRB.fa`. Therefore, to get a [detailed TSV file](/manual/#output-files) with per-read clonotype information, **IMSEQ** can be invoked as follows:

~~~Plaintext
$ imseq -ref Homo.Sapiens.TRB.fa -o result.tsv -r examples/data/example_sim.fa
~~~

**Note the additional use of `-r`!** The default behaviour of **IMSEQ** is to consider the V read file to be forward sequence and the VDJ read file to be complementary reverse. While this reflects a common scenario in real world experiments, the provided simulated sequences are not complementary reverse.

The specified output file `result.tsv` should look like this:

~~~Plaintext
seqId    cdrBegin  cdrEnd  leftMatches  leftErrPos  leftMatchLen  rightMatches  rightErrPos  rightMatchLen  cdrNucSeq                                cdrAASeq
READ-1   273       300     V16                      276           J2-7                       31             TGTGCCAGCACCTACGAGCAGTACTTC              CASTYEQYF
READ-10  270       303     V6-2V6-3                 273           J2-4                       31             TGTGCCAGCATTTCCAAAGACATTCAGTACTTC        CASISKDIQYF
READ-11  270       297     V3-1                     273           J2-1                       31             TGTGCCTGCTACAATGAGCAGTTCTTC              CACYNEQFF
READ-12  273       303     V7-4                     276           J2-1                       31             TGTGCCGGCTGCTTAGATGAGCAGTTCTTC           CAGCLDEQFF
READ-13  273       312     V23-1                    276           J2-6                       31             TGCGCCAGCAGGCAATCTGAGGCCAACGTCCTGACTTTC  CASRQSEANVLTF
READ-14  257       284     V11-3                    260           J1-4                       31             TGTGCCACCAGTGAAAAACTGTTTTTT              CATSEKLFF
READ-15  270       297     V6-6                     273           J2-4                       31             TGTGCCGGCAGAAACATTCAGTACTTC              CAGRNIQYF
READ-16  270       300     V15                      273           J2-7                       31             TGTGCCACCATCAGCGACGAGCAGTACGTC           CATISDEQYV
READ-17  273       309     V7-2                     276           J1-5                       31             TGTGCCAGCAGCTTTGGGCATCAGCCCCAGCATTTT     CASSFGHQPQHF
...
~~~
