---
layout: default
permalink: manual/
---
# Manual

* This will become a table of contents (this text will be scraped).
{:toc}


## Introduction

**IMSEQ** is a tool for the analysis of T- and B-cell receptor chain sequences. It can be used to analyse either single-read data, where the reads cover the V-CDR3-J region sufficiently for an identification, or paired-end data where one read covers the V-region and one read covers the J- and CDR3-region. The latter read has to cover only a small fraction of the V-segment, sufficient for the localization of the Cys-104 motif (TODO- Correct spelling in **IMSEQ** app).  This manual will begin with a listing of the command line arguments and their intended use, and then preset some vignettes that demonstrate the intended use of **IMSEQ**. The simplest **IMSEQ** command is used to show a complete list of options on the com- mand line:

    $ imseq --help

This command will show the expected syntax of each of the command line arguments, e.g.:

    -oa, --out-amino STR
          Output file path for translated clonotypes.

For each argument, there is a short form (here: -oa) which is predeced by a single dash, and an equivalent long form (--out-amino) that is preced by two dashes. Either form can be used interchangably. For each argument, the required parameter type is given:

| STR  | A string, i.e., and combination of letters, numerals, or the underscore character (\_) |
| NUM  | A valid number |
| FILE | Absolute or relative path to a file that already exists on the file system. |

<p></p>
If no parameter type is specified, the argument does not accept any parameters and is used as a flag to enable or disable a feature.

## Basic usage and input files

For single-end V(D)J-reads, **IMSEQ** is called as follows:

    $ imseq -ref <segment sequences> {-o,-oa,-on} <output file> <VDJ-reads>

For paired-end data, the program call takes two input files:

    $ imseq -ref <segment sequences> {-o,-oa,-on} <output file> <V-reads> <VDJ-reads>

The segment sequences have to be provided in FASTA format as specified in the [IMSEQ FASTA ID Specification](../fastaFormat). At least one output file has to be provided, using the -o, -oa or -on option. For multiple output files, each option has to be provided with its own output file name. The input sequences / reads can be provided either in FASTA or FASTQ format, while all quality related features only work when FASTQ files were provided. The input files can be GZIP compressed.

## Output files

**IMSEQ** supports three output files to write the clonotyping and repertoire generation results:

 * ```-o, --out``` - _Detailed per read clonotyping information_<br/>This file contains one line for every read that was successfully analysed, providing the following information for every input sequences that was successfully analyzed:
   1. _seqId_<br/>The FASTA / FASTQ ID of the input sequence
   1. _cdrBegin_<br/>The position of the first base of the Cys-triplet in the read, 0-based
   1. _cdrEnd_<br/>The position of the last base of the Phe-triplet in the read, 0-based
   1. _leftMatches_<br/>The identified V segment(s)
   1. _leftErrPos_<br/>The read to segment alignment mismatch positions relative to the last base before the Cys-triplet
   1. _leftMatchLen_<br/>The length of the V segment alignment until the last base of the segment core fragment
   1. _rightMatches_<br/>The identified J segment(s)
   1. _rightErrPos_<br/>The read to segment alignment mismatch positions relative to the first base after the Phe-triplet
   1. _rightMatchLen_<br/>The length of the V segment alignment from the first base of the segment core fragment on
   1. _cdrNucSeq_<br/>The nucleotide sequence of the CDR3 region
   1. _cdrAASeq_<br/>The aminoacid sequence of the CDR3 region
 * ```-on, --out-nuc``` - _Nucleotide based clonotype counts_<br/>Counts of clonotypes, where a clonotype is defined by its V segment, J segment and the nucleotide sequence of the CDR3 region. The counts are generated after the correction of PCR and sequencing errors.
 * ```-oa, --out-amino``` - _Aminoacid based clonotype counts_<br/>Potentially corrected counts of clonotypes, where a clonotype is defined by its V segment, J segment and the aminoacid sequence of the CDR3 region. The counts are generated after the correction of PCR and sequencing errors.

Furthermore, the user can choose to specify one of the following options:

 * ```-al, --with-alleles``` - Split clonotypes into different V/J alleles for ```-oa``` and ```-on```.
 * ```-rlog, --reject-log``` - _A list of IDs of rejected reads_<br/>For every rejected read, the reason for the rejection is given as one of the following:<br/>

    | ```AVERAGE_QUAL_FAIL``` | The average base quality score of the read was too low (see ```-mq```) |
    | ```MOTIF_AMBIGUOUS``` | Equally well matching segments imply different Cys/Phe triplet position |
    | ```NONSENSE_IN_CDR3``` | CDR3 region contains STOP codon |
    | ```OUT_OF_READING_FRAME``` | Cys and Phe are out of frame |
    | ```SEGMENT_MATCH_FAILED``` | No V or no J segment could be identified |

## Read preprocessing

 * ```-r, --reverse``` - _Change default reverse complementation behaviour_<br/>By default, **IMSEQ** reverse complements the V(D)J-reads. If this option is specified, the V-reads are reverse complemented instead.
 * ```-tr, --truncate-reads``` - _Truncate reads to specified length_<br/>Removes the specified number of bases from the end of the input reads before processing them further.

## V/J segment alignment

 * ```-ev, --v-err-rate NUM``` - _Maximum error rate for the V segment to V(D)J read alignment._<br/>The maximum error rate allowed for matching a V segment against the V(D)J-read. Default: 0.05.
 * ```-ej, --j-err-rate NUM``` - _Maximum error rate for the J segment to V(D)J read alignment._<br/>The maximum error rate allowed for matching a J segment against the V(D)J-read. Default: 0.15.

## V/J segment alignment (paired-end)

  * ```-pve, --paired-v-error``` - _Maximum error rate for the V segment to V read alignment._<br/>The maximum error rate allowed for matching a V segment against the V-read. Default: Use value from ```-ev```.

## V/J segment alignment (Expert settings)

  * ```-jcl, --j-core-length``` - _Length of the J core fragment._
  * ```-jco, --j-core-offset``` - _Offset of the V core fragment._
  * ```-vcl, --v-core-length``` - _Length of the V core fragment._
  * ```-vco, --v-core-offset``` - _Offset of the V core fragment._
  * ```-vce, --v-core-errors``` - _Maximum number of errors when matching the V core fragments._
  * ```-jce, --j-core-errors``` - _Maximum number of errors when matching the J core fragments._

## Quality filtering

  * ```-mq, --min-qual``` - _Minimum average read phred score._<br/>The minimum average base quality score required for a read to be analysed. Default: 10.
  * ```-mcq, --min-clust-qual``` - _Minimum average cluster phred score._<br/>The minimum average read phred score across an entire clonotype cluster. This filter is applied after clustering and can therefore be used to remove low quality clonotype cluster that couldn't be corrected. Default: 30.

## Postprocessing / Clustering

  * ```-ma, --merge-ambiguous-seg``` - _Enable segment ambiguity clustering_<br/>Segment ambiguity clustering reduces the number of clonotypes with ambiguous segment assignments by merging them with clonotypes with less or no ambiguity in the segment assignments that match with respect to the CDR3 sequence.
  * ```-qc, --qual-clust``` - _Enable quality score based clustering._
  * ```-sc, --simple-clust``` - _Enable simple distance-based clustering_
  * ```-qcme, --max-err-hq``` - _Maximum edit-distance for two clusters to be clustered without low quality correlation_
  * ```-qcsd, --min-sd-diff``` - _How many standard deviations must an error positions quality value be below the mean to be considered for_
  * ```-scme, --max-err-lq``` - _Maximum edit-distance for two clusters to be potentially clustered without low quality correlation_
  * ```-mcr, --max-clust-ratio``` - _Maximum abundance ratio for two clonotypes to be clustered_<br/>Sets the maximum ratio A/B, where A is the minor clonotypes count and B is the major clonotypes count. Default: no restriction.

## Performance

  * ```-j, --jobs``` - _Number of parallel jobs (threads)_<br/>Default: 1.

## Other options

 * ```-pa, --print-alignments``` - _Print the V/J alignments for each read. Implies -j 1._
