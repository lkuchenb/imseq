// ============================================================================
// IMSEQ - An immunogenetic sequence analysis tool
// (C) Charite, Universitaetsmedizin Berlin
// Author: Leon Kuchenbecker
// ============================================================================
// 
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published by
// the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
// 
// You should have received a copy of the GNU General Public License along with
// this program; if not, write to the Free Software Foundation, Inc., 51
// Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// ============================================================================

// ============================================================================
// HEADER FILE DESCRIPTION
// ============================================================================
// Reading FASTA / FASTQ files
// ============================================================================

#ifndef IMSEQ_FASTQ_IO_TYPES_H
#define IMSEQ_FASTQ_IO_TYPES_H

#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "sequence_data_types.h"

/********************************************************************************
 * STRUCTS AND CLASSES
 ********************************************************************************/

/**
 * A FastqRecord consists either of one or of two sequences (single-end,
 * paired-end), an ID (only one ID is supported) and a barcode sequence.
 */
template<typename T>
struct FastqRecord {};

/**
 * @special Paired end implementation
 */
template<>
struct FastqRecord<PairedEnd> {
    typedef String<Dna5Q> TSequence;
    CharString  id;
    TSequence   fwSeq, revSeq, bcSeq;
    FastqRecord<PairedEnd>() {}
    FastqRecord<PairedEnd>(CharString const & _id, TSequence const & _fwSeq, TSequence const & _revSeq, TSequence const & _bcSeq) :
        id(_id), fwSeq(_fwSeq), revSeq(_revSeq), bcSeq(_bcSeq) {}
};

/**
 * @special Single end implementation
 */
template<>
struct FastqRecord<SingleEnd> {
    typedef String<Dna5Q> TSequence;
    CharString  id;
    TSequence   seq, bcSeq;

    FastqRecord<SingleEnd>() {}
    FastqRecord<SingleEnd>(CharString const & _id, TSequence const & _seq, TSequence const & _bcSeq) :
        id(_id), seq(_seq), bcSeq(_bcSeq) {}
    FastqRecord<SingleEnd> (FastqRecord<PairedEnd> const & pairRec)
    {
	SEQAN_CHECK(empty(pairRec.fwSeq), "PLEASE REPORT THIS ERROR");
	seq = pairRec.revSeq;
	id  = pairRec.id;
    }
};

inline void openOrExit(SeqFileIn & stream, std::string const & path)
{
    if (!open(stream, path.c_str())) {
        std::cerr << "[ERROR] Cannot open '" << path << "'" << std::endl;
        std::exit(1);
    }
}

/**
 * SeqInputStreams can hold either one or two paths to sequence files and the
 * corresponding SeqFileIn objects. Upon construction with given path(s) it
 * opens the streams.
 */
template<typename TSequencingType>
struct SeqInputStreams {};

/**
 * @special Single end implementation
 */
template<>
struct SeqInputStreams<SingleEnd> {
    std::string path;
    SeqFileIn stream;
    uint64_t totalInBytes;
    SeqInputStreams<SingleEnd>(std::string path_) : path(path_), totalInBytes(0) {
        openOrExit(stream, path.c_str());
    }
    SeqInputStreams<SingleEnd>(std::string path_, uint64_t totalInBytes_) : path(path_), totalInBytes(totalInBytes_) {
        openOrExit(stream, path.c_str());
    }
};

/**
 * @special Paired end implementation
 */
template<>
struct SeqInputStreams<PairedEnd> {
    std::string fwPath, revPath;
    SeqFileIn fwStream, revStream;
    uint64_t totalInBytes;
    SeqInputStreams<PairedEnd>(std::string fwPath_, std::string revPath_) : fwPath(fwPath_), revPath(revPath_), totalInBytes(0) {
        openOrExit(fwStream,fwPath.c_str());
        openOrExit(revStream,revPath.c_str());
    }
    SeqInputStreams<PairedEnd>(std::string fwPath_, std::string revPath_, uint64_t totalInBytes_) : fwPath(fwPath_), revPath(revPath_), totalInBytes(totalInBytes_) {
        openOrExit(fwStream,fwPath.c_str());
        openOrExit(revStream,revPath.c_str());
    }
};

#endif
