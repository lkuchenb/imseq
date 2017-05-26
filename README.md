# IMSEQ - An immunogenetic sequence analysis tool

[![Build Status](https://travis-ci.org/lkuchenb/imseq.svg?branch=master)](https://travis-ci.org/lkuchenb/imseq) [![GitHub release](https://img.shields.io/github/release/lkuchenb/imseq.svg)](https://github.com/lkuchenb/imseq/releases/latest)

This is the source code repository for **IMSEQ**. For information about **IMSEQ** and release downloads visit www.imtools.org!

## Building **IMSEQ**

These instructions apply if you want to build **IMSEQ** from the git repository.

### Requirements

 - A C++14 capable compiler
 - ZLIB

### Instructions

Clone the repository

    $ git clone https://github.com/lkuchenb/imseq.git
    $ cd imseq

Configure a build directory

    $ mkdir build
    $ cd build
    $ cmake ..

Build **IMSEQ**

    $ make imseq

### Troubleshooting

 - Cmake will attempt to download the SeqAn Library. If this fails or you want to use your own copy of SeqAn, invoke cmake with `-DSEQAN_ROOT=/path/to/seqan`.
 - If you want to use a compiler different from your systems default compiler, invoke cmake with `-DCMAKE_CXX_COMPILER=/path/to/c++`
