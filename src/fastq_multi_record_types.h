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
// ============================================================================

#ifndef IMSEQ_FASTQ_MULTI_RECORD_TYPES_H
#define IMSEQ_FASTQ_MULTI_RECORD_TYPES_H

#include <map>
#include <unordered_map>
#include <limits>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "sequence_data_types.h"

using namespace seqan;

/*-------------------------------------------------------------------------------
 - FastqMultiRecord
 -------------------------------------------------------------------------------*/

/**
 * A FastqMultiRecord holds reads that share the same sequence
 */
template<typename T>
struct FastqMultiRecord {};

/**
 * @special Single end implementation
 */
template<>
struct FastqMultiRecord<SingleEnd> {
    typedef String<Dna5> TSequence;
    typedef String<float> TQualities;
    typedef std::set<CharString> TIds;

    TSequence   seq, bcSeq;
    TQualities  qualities;
    TIds        ids;
};

/**
 * @special Paired end implementation
 */
template<>
struct FastqMultiRecord<PairedEnd> {
    typedef String<Dna5> TSequence;
    typedef String<float> TQualities;
    typedef std::set<CharString> TIds;

    TSequence   fwSeq, revSeq, bcSeq;
    TQualities  fwQualities, revQualities;
    TIds        ids;
};

/*-------------------------------------------------------------------------------
 - FastqMultiRecordCollection
 -------------------------------------------------------------------------------*/

struct DnaStringHash
{
    static std::hash<std::string> const hash_fun;
    size_t operator()(String<Dna5> const & str) const;
};

template<typename T>
struct FastqMultiRecordCollection {};

template<>
struct FastqMultiRecordCollection<SingleEnd> {
    static uint64_t const NO_MATCH = std::numeric_limits<uint64_t>::max();

//    typedef std::unordered_map<FastqMultiRecord<SingleEnd>::TSequence, unsigned, DnaStringHash> TSeqMap;
//    typedef std::unordered_map<FastqMultiRecord<SingleEnd>::TSequence, TSeqMap, DnaStringHash> TBcMap;

    typedef std::map<FastqMultiRecord<SingleEnd>::TSequence, unsigned> TSeqMap;
    typedef std::map<FastqMultiRecord<SingleEnd>::TSequence, TSeqMap> TBcMap;

    String<FastqMultiRecord<SingleEnd> > multiRecords;
    TBcMap bcMap;
};

template<>
struct FastqMultiRecordCollection<PairedEnd> {
    static uint64_t const NO_MATCH = std::numeric_limits<uint64_t>::max();

//    typedef std::unordered_map<FastqMultiRecord<PairedEnd>::TSequence, unsigned, DnaStringHash> TRevSeqMap;
//    typedef std::unordered_map<FastqMultiRecord<PairedEnd>::TSequence, TRevSeqMap, DnaStringHash> TFwSeqMap;
//    typedef std::unordered_map<FastqMultiRecord<PairedEnd>::TSequence, TFwSeqMap, DnaStringHash> TBcMap;

    typedef std::map<FastqMultiRecord<PairedEnd>::TSequence, unsigned> TRevSeqMap;
    typedef std::map<FastqMultiRecord<PairedEnd>::TSequence, TRevSeqMap> TFwSeqMap;
    typedef std::map<FastqMultiRecord<PairedEnd>::TSequence, TFwSeqMap> TBcMap;

    String<FastqMultiRecord<PairedEnd> > multiRecords;
    TBcMap bcMap;
};

/*-------------------------------------------------------------------------------
 - BarcodeStats
 -------------------------------------------------------------------------------*/

struct BarcodeStats {
    typedef StringSet<String<Dna5> > TBcSeqs;
    TBcSeqs bcSeqs;
    String<uint64_t> nReads;
    String<uint64_t> nUniqueReads;
    uint64_t nTotalUniqueReads, nTotalReads;
    BarcodeStats() : nTotalUniqueReads(0), nTotalReads(0) {}
};

#endif
