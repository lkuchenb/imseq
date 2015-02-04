---
layout: default
date:   2015-02-04 15:53:59
permalink: fastaFormat/
---
# FASTA ID Specification

The V and J segment reference sequences must be provided in FASTA format. The FASTA sequence IDs must contain the following information, separated by ```|``` characters:

 * The gene, i.e. TRB, IGH, etc
 * The segment type, either ```V``` or ```J```
 * The segment ID
 * The allele ID
 * The position of the first base of the Cys (or Phe) triplet (counting from 0)

Both the V and J segments have to be provided in the same file. Such a file might look like this:

    ...
    >TRB|V|9|02|270
    gattctggagtcacacaaaccccaaagcacctgatcacagcaactggacagcgagtgacg
    ctgagatgctcccctaggtctggagacctctctgtgtactggtaccaacagagcctggac
    cagggcctccagttcctcattcactattataatggagaagagagagcaaaaggaaacatt
    cttgaacgattctccgcacaacagttccctgacttgcactctgaactaaacctgagctct
    ctggagctgggggactcagctttgtatttctgtgccagcagcgtag
    >TRB|V|9|03|270
    gattctggagtcacacaaaccccaaagcacctgatcacagcaactggacagcgagtgacg
    ctgagatgctcccctaggtctggagacctctctgtgtactggtaccaacagagcctggac
    cagggcctccagttcctcattcaatattataatggagaagagagagcaaaaggaaacatt
    cttgaacgattctccgcacaacagttccctgacttgcactctgaactaaacctgagctct
    ctggagctgggggactcagctttgtatttctgtgccagcagc
    >TRB|J|1-1|01|17
    tgaacactgaagctttctttggacaaggcaccagactcacagttgtag
    >TRB|J|1-2|01|17
    ctaactatggctacaccttcggttcggggaccaggttaaccgttgtag
    ...
