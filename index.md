---
layout: default
---

IMSEQ is a **fast, PCR and sequencing error aware** tool to analyze high throughput data from recombined **T-cell receptor** or **immunoglobolin** gene **sequencing** experiments. It derives immune **repertoires** from sequencing data in FASTA / FASTQ format.

![IMSEQ Workflow](./images/imseq-flow.png)

## Getting IMSEQ

You can obtain binaries and sources for the latest release, **IMSEQ 1.0**, by choosing from one of the following links:

 * [**IMSEQ 1.0 for Linux (64 bit)**](./binaries/imseq_1.0-linux64.tgz)<br/>Statically linked binary supporting most 64 bit Linux systems
 * [**IMSEQ 1.0 for OS X (64 bit)**](./binaries/imseq_1.0-mac64.tgz)<br/>Requires OS X 10.9 or newer
 * [**IMSEQ 1.0 Sourcecode**](./binaries/seqan-imseq_1.0-source.tgz)<br/>A source tarball of IMSEQ 1.0. For the sources and build instructions of the current development version check out the [GitHub repository](https://github.com/lkuchenb/imseq).
 
## Using IMSEQ

**IMSEQ** requires at least two input files and the specification of at least one output file. The two input files are specified using the switches

    $ imseq -ref segment-reference.fa -o output-file.tsv input-file.fastq.gz

<p></p>
 * _segment-reference.fa_<br/>A file containing the ***V and J segment sequences*** for the gene and species that is analyzed. The file must be in FASTA format and the sequence IDs must follow the [IMSEQ FASTA ID Specification](./fastaFormat). Gene segment reference files for the human T-cell receptor alpha and beta chain genes are provided with **IMSEQ**.
 * _input-file.fastq.gz_<br/>One (single-end sequencing) or two (split paired-end sequencing) **FASTA or FASTQ files** with the **input reads**. If two files are specified, the first file is considered to be V-read and the second file the V(D)J-read, if only one file is specified it has to contain V(D)J-reads. V-reads are parsed as forward sequence, V(D)J-reads are parsed as reverse complementary sequence. This behavior can be inversed using the ```-r``` switch.
 * _output-file.tsv_<br/>An output file where the detailed per-read analysis results are written in TSV format. See the [output file specifications](./manual#output-files) for more details.

## Further Reading

**IMSEQ** supports many more options that can adjust its behaviour to many experimental conditions. All options are documented in the [Manual](./manual). Furthermore, the [Tutorial](./tutorial) will guide you through some basic analysis of the example files that come with **IMSEQ**.
