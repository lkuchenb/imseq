---
layout: default
---

IMSEQ is a **fast, PCR and sequencing error aware** tool to analyze high throughput data from recombined **T-cell receptor** or **immunoglobolin** gene **sequencing** experiments. It derives immune **repertoires** from sequencing data in FASTA / FASTQ format.

![IMSEQ Workflow](./images/imseq-flow.png)

<div style="border:1pt solid black; padding:5px; box-shadow:5px 5px grey; margin-bottom:15px" markdown="block">
Please cite the following publication when you use **IMSEQ**:

Kuchenbecker L, Nienen M, Hecht J, Neumann AU, Babel N, Reinert K, Robinson PN. *IMSEQ - a fast and error aware approach to immunogenetic sequence analysis*. Bioinformatics. 2015;31(18):2963–71.
[\[PubMed\]](http://www.ncbi.nlm.nih.gov/pubmed/25987567) [\[Journal\]](http://bioinformatics.oxfordjournals.org/content/31/18/2963) [\[BibTeX\]](./bibcite)
</div>

## Getting IMSEQ

**IMSEQ** is freely available and its source code is made available under the [GPLv2](http://www.gnu.org/licenses/gpl-2.0.html) license. You can obtain binaries and sources for the latest release, **IMSEQ 1.0.3**, by choosing from one of the following links:

 * [**IMSEQ 1.0.3 for Linux (64 bit)**](https://github.com/lkuchenb/imseq/releases/download/v1.0.3/imseq_1.0.3-linux64.tgz)<br/>Statically linked binary supporting most 64 bit Linux systems
 * [**IMSEQ 1.0.3 for OS X (64 bit)**](https://github.com/lkuchenb/imseq/releases/download/v1.0.3/imseq_1.0.3-mac64.tgz)<br/>Requires OS X 10.9 or newer
 * [**IMSEQ 1.0.3 Sourcecode**](https://github.com/lkuchenb/imseq/releases/download/v1.0.3/seqan-imseq_1.0.3-source.tgz)<br/>A source tarball of IMSEQ 1.0.3. 

Changes can be reviewed in the [**Changelog**](https://github.com/lkuchenb/imseq/releases/).

For the sources and build instructions of the current **development version** check out the [GitHub repository](https://github.com/lkuchenb/imseq).
 
## Using IMSEQ

**IMSEQ** requires at least two input files and the specification of at least one output file. The two input files are specified using the switches

    $ imseq -ref segment-reference.fa -o output-file.tsv input-file.fastq.gz

<p></p>
 * _segment-reference.fa_<br/>A file containing the ***V and J segment sequences*** for the gene and species that is analyzed. The file must be in FASTA format and the sequence IDs must follow the [IMSEQ FASTA ID Specification](./fastaFormat). Gene segment reference files for the human T-cell receptor alpha and beta chain as well as immunoglobolin heavy, light-kappa and light-lambda genes are provided with **IMSEQ**.
 * _input-file.fastq.gz_<br/>One (single-end sequencing) or two (split paired-end sequencing) **FASTA or FASTQ files** with the **input reads**. If two files are specified, the first file is considered to be V-read and the second file the V(D)J-read, if only one file is specified it has to contain V(D)J-reads. V-reads are parsed as forward sequence, V(D)J-reads are parsed as reverse complementary sequence. This behavior can be inversed using the ```-r``` switch.
 * _output-file.tsv_<br/>An output file where the detailed per-read analysis results are written in TSV format. See the [output file specifications](./manual#output-files) for more details.

## Further Reading

**IMSEQ** supports many more options that can adjust its behaviour to many experimental conditions. All options are documented in the [Manual](./manual). Furthermore, the [Tutorial](./tutorial) will guide you through some basic analysis of the example files that come with **IMSEQ**.
