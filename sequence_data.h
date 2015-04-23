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

#ifndef CDR3FINDER_SEQUENCE_DATA_H
#define CDR3FINDER_SEQUENCE_DATA_H

#include <seqan/basic.h>
#include "overlap_specs.h"

/********************************************************************************
 * TAGS
 ********************************************************************************/

struct SingleEnd_;
typedef Tag<SingleEnd_> SingleEnd;

struct PairedEnd_;
typedef Tag<PairedEnd_> PairedEnd;

/********************************************************************************
 * STRUCTS AND CLASSES
 ********************************************************************************/

typedef String<Dna5Q>           TQueryDataSequence;

template<typename T>
struct QueryData {
    typedef String<Dna5Q>	TSequence;
};

template<>
struct QueryData<SingleEnd> : QueryData<void> {
    StringSet<TQueryDataSequence>       seqs;
    StringSet<CharString>               ids;
};

template<>
struct QueryData<PairedEnd> : QueryData<void> {
    StringSet<TQueryDataSequence>       fwSeqs, revSeqs;
    StringSet<CharString>               ids;
};

template<typename TSequencingSpec>
struct QueryDataCollection {};

template<>
struct QueryDataCollection<SingleEnd>
{
    QueryData<SingleEnd> queryData;
};

template<>
struct QueryDataCollection<PairedEnd>
{
    QueryData<PairedEnd> pairedQueryData;
    QueryData<SingleEnd> singleQueryData;
};

template<typename T>
struct FastqRecord {};

template<>
struct FastqRecord<PairedEnd> {
    typedef String<Dna5Q> TSequence;
    TSequence   fwSeq, revSeq;
    CharString  id;
};

template<>
struct FastqRecord<SingleEnd> {
    typedef String<Dna5Q> TSequence;
    TSequence   seq;
    CharString  id;

    FastqRecord<SingleEnd> () {}
    FastqRecord<SingleEnd> (FastqRecord<PairedEnd> const & pairRec)
    {
	SEQAN_CHECK(empty(pairRec.fwSeq), "PLEASE REPORT THIS ERROR");
	seq = pairRec.revSeq;
	id  = pairRec.id;
    }
};

/********************************************************************************
 * FUNCTIONS
 ********************************************************************************/

inline void printRecord(FastqRecord<PairedEnd> const & rec) {
        std::cerr << "FORWARD\t" << rec.fwSeq
            << "\nREVERSE\t" << rec.revSeq << '\n';
}

inline void printRecord(FastqRecord<SingleEnd> const & rec) {
        std::cerr << "READ\t" << rec.seq << '\n';
}

inline unsigned nRecords(QueryData<SingleEnd> const & qData)
{
    return length(qData.seqs);
}

inline unsigned nRecords(QueryData<PairedEnd> const & qData)
{
    return length(qData.fwSeqs);
}

inline unsigned nRecords(QueryDataCollection<SingleEnd> const & qDatCol) {
    return nRecords(qDatCol.queryData);
}

inline unsigned nRecords(QueryDataCollection<PairedEnd> const & qDatCol) {
    return nRecords(qDatCol.singleQueryData) + nRecords(qDatCol.pairedQueryData);
}

inline void appendRecord(QueryData<SingleEnd> & target, QueryData<SingleEnd> const & source, unsigned pos) {
    appendValue(target.seqs, source.seqs[pos]);
    appendValue(target.ids, source.ids[pos]);
}

inline void appendRecord(QueryData<PairedEnd> & target, QueryData<PairedEnd> const & source, unsigned pos) {
    appendValue(target.fwSeqs, source.fwSeqs[pos]);
    appendValue(target.revSeqs, source.revSeqs[pos]);
    appendValue(target.ids, source.ids[pos]);
}

// ============================================================================
// Getter for V Read Sequence
// ============================================================================

inline StringSet<TQueryDataSequence> & getVReadSequences(QueryData<SingleEnd> & qData)
{
    return qData.seqs;
}

inline StringSet<TQueryDataSequence> & getVReadSequences(QueryData<PairedEnd> & qData)
{
    return qData.fwSeqs;
}

inline StringSet<TQueryDataSequence> const & getVReadSequences(QueryData<SingleEnd> const & qData)
{
    return qData.seqs;
}

inline StringSet<TQueryDataSequence> const & getVReadSequences(QueryData<PairedEnd> const & qData)
{
    return qData.fwSeqs;
}

// ============================================================================
// Getter for VDJ Read Sequence
// ============================================================================

inline StringSet<TQueryDataSequence> & getVDJReadSequences(QueryData<SingleEnd> & qData)
{
    return qData.seqs;
}

inline StringSet<TQueryDataSequence> & getVDJReadSequences(QueryData<PairedEnd> & qData)
{
    return qData.revSeqs;
}

inline StringSet<TQueryDataSequence> const & getVDJReadSequences(QueryData<SingleEnd> const & qData)
{
    return qData.seqs;
}

inline StringSet<TQueryDataSequence> const & getVDJReadSequences(QueryData<PairedEnd> const & qData)
{
    return qData.revSeqs;
}

// ============================================================================
// Getter for V Read IDs
// ============================================================================

template<typename TSequencingSpec>
StringSet<CharString> & getReadIds(QueryData<TSequencingSpec> & qData)
{
    return qData.ids;
}

template<typename TSequencingSpec>
StringSet<CharString> const & getReadIds(QueryData<TSequencingSpec> const & qData)
{
    return qData.ids;
}

// ============================================================================
// Overlap Tag based functions
// ============================================================================

template<typename TQueryData>
StringSet<TQueryDataSequence> & getReadSequences(TQueryData & qData, LeftOverlap const) 
{
    return getVReadSequences(qData);
}

template<typename TQueryData>
StringSet<TQueryDataSequence> & getReadSequences(TQueryData & qData, RightOverlap const) 
{
    return getVDJReadSequences(qData);
}

template<typename TQueryData>
StringSet<TQueryDataSequence> const & getReadSequences(TQueryData const & qData, LeftOverlap const) 
{
    return getVReadSequences(qData);
}

template<typename TQueryData>
StringSet<TQueryDataSequence> const & getReadSequences(TQueryData const & qData, RightOverlap const) 
{
    return getVDJReadSequences(qData);
}

template<typename TQueryData>
StringSet<CharString> & getReadIds(TQueryData & qData, LeftOverlap const) 
{
    return getVReadIds(qData);
}

template<typename TQueryData>
StringSet<CharString> & getReadIds(TQueryData & qData, RightOverlap const) 
{
    return getVDJReadIds(qData);
}

template<typename TQueryData>
StringSet<CharString> const & getReadIds(TQueryData const & qData, LeftOverlap const) 
{
    return getVReadIds(qData);
}

template<typename TQueryData>
StringSet<CharString> const & getReadIds(TQueryData const & qData, RightOverlap const) 
{
    return getVDJReadIds(qData);
}

#endif
