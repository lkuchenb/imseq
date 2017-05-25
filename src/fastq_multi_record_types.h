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

#include <unordered_map>
#include <list>
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
    typedef String<double> TQualities;
    typedef std::set<CharString> TIds;

    TSequence   seq, bcSeq;
    std::set<TSequence> bcSeqHistory;
    TQualities  qualities;
    TIds        ids;
};

/**
 * @special Paired end implementation
 */
template<>
struct FastqMultiRecord<PairedEnd> {
    typedef String<Dna5> TSequence;
    typedef String<double> TQualities;
    typedef std::set<CharString> TIds;

    TSequence   fwSeq, revSeq, bcSeq;
    std::set<TSequence> bcSeqHistory;
    TQualities  fwQualities, revQualities;
    TIds        ids;
};

/*-------------------------------------------------------------------------------
 - FastqMultiRecordCollection
 -------------------------------------------------------------------------------*/

namespace std {

  template<>
  struct hash<seqan::String<Dna5> >
  {
    std::size_t operator()(const seqan::String<Dna5> & k) const
    {
        size_t res = 12345678910;
        for (seqan::Iterator<seqan::String<Dna5> const,Rooted>::Type it = begin(k); !atEnd(it); goNext(it))
        {
            res = res * 101 + static_cast<int>(*it) * 60;
        }
        return res;
    }
  };

}

template<typename T>
struct FastqMultiRecordCollection {};

template<>
struct FastqMultiRecordCollection<SingleEnd> {
    typedef FastqMultiRecord<SingleEnd>         TMRec;
    typedef TMRec::TSequence                    TSequence;
    typedef std::list<TMRec>                    TRecList;
    typedef std::list<TMRec>::const_iterator    TRecListIt;

    static uint64_t const NO_MATCH = std::numeric_limits<uint64_t>::max();

    typedef std::unordered_map<TMRec::TSequence, TMRec*> TSeqMap;
    typedef std::unordered_map<TMRec::TSequence, TSeqMap> TBcMap;

    TRecList multiRecords;
    TBcMap bcMap;
};

template<>
struct FastqMultiRecordCollection<PairedEnd> {
    typedef FastqMultiRecord<PairedEnd>         TMRec;
    typedef TMRec::TSequence                    TSequence;
    typedef std::list<TMRec>                    TRecList;
    typedef std::list<TMRec>::const_iterator    TRecListIt;

    static uint64_t const NO_MATCH = std::numeric_limits<uint64_t>::max();

    typedef std::unordered_map<TMRec::TSequence, TMRec*> TRevSeqMap;
    typedef std::unordered_map<TMRec::TSequence, TRevSeqMap> TFwSeqMap;
    typedef std::unordered_map<TMRec::TSequence, TFwSeqMap> TBcMap;

    TRecList multiRecords;
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
