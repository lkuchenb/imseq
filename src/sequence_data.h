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
#include <seqan/sequence.h>

#include "overlap_specs.h"
#include "fastq_io_types.h"
#include "sequence_data_types.h"


/********************************************************************************
 * FUNCTIONS
 ********************************************************************************/

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

//inline void appendRecord(QueryData<SingleEnd> & target, QueryData<SingleEnd> const & source, unsigned pos) {
//    appendValue(target.seqs, source.seqs[pos]);
//    appendValue(target.bcSeqs, source.bcSeqs[pos]);
//    appendValue(target.ids, source.ids[pos]);
//}
//
//inline void appendRecord(QueryData<PairedEnd> & target, QueryData<PairedEnd> const & source, unsigned pos) {
//    appendValue(target.fwSeqs, source.fwSeqs[pos]);
//    appendValue(target.revSeqs, source.revSeqs[pos]);
//    appendValue(target.bcSeqs, source.bcSeqs[pos]);
//    appendValue(target.ids, source.ids[pos]);
//}
//
inline void clear(QueryData<SingleEnd> & qData) {
    clear(qData.seqs);
}

inline void clear(QueryData<PairedEnd> & qData) {
    clear(qData.fwSeqs);
    clear(qData.revSeqs);
}

//inline void addRecord(QueryData<PairedEnd> & qData, FastqRecord<PairedEnd> const & fqRecord)
//{
//    appendValue(qData.fwSeqs, fqRecord.fwSeq);
//    appendValue(qData.revSeqs, fqRecord.revSeq);
//    appendValue(qData.bcSeqs, fqRecord.bcSeq);
//    appendValue(qData.ids, fqRecord.id);
//}
//
//inline void addRecord(QueryData<SingleEnd> & qData, FastqRecord<SingleEnd> const & fqRecord)
//{
//    appendValue(qData.seqs, fqRecord.seq);
//    appendValue(qData.bcSeqs, fqRecord.bcSeq);
//    appendValue(qData.ids, fqRecord.id);
//}

//inline void addRecord(QueryDataCollection<SingleEnd> & qDatCol, FastqRecord<SingleEnd> const & rec)
//{
//    addRecord(qDatCol.queryData, rec);
//}
//
//inline void addRecord(QueryDataCollection<PairedEnd> & qDatCol, FastqRecord<PairedEnd> const & rec)
//{
//    if (empty(rec.fwSeq)) 
//        addRecord(qDatCol.singleQueryData, rec);
//    else
//        addRecord(qDatCol.pairedQueryData, rec);
//}

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

inline StringSet<String<double> > & getVDJAvgQuals(QueryData<SingleEnd> & qData)
{
    return qData.avgQVals;
}

inline StringSet<String<double> > & getVDJAvgQuals(QueryData<PairedEnd> & qData)
{
    return qData.revAvgQVals;
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

//template<typename TSequencingSpec>
//StringSet<CharString> & getReadIds(QueryData<TSequencingSpec> & qData)
//{
//    return qData.ids;
//}
//
//template<typename TSequencingSpec>
//StringSet<CharString> const & getReadIds(QueryData<TSequencingSpec> const & qData)
//{
//    return qData.ids;
//}
//
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

//template<typename TQueryData>
//StringSet<CharString> & getReadIds(TQueryData & qData, LeftOverlap const) 
//{
//    return getVReadIds(qData);
//}
//
//template<typename TQueryData>
//StringSet<CharString> & getReadIds(TQueryData & qData, RightOverlap const) 
//{
//    return getVDJReadIds(qData);
//}
//
//template<typename TQueryData>
//StringSet<CharString> const & getReadIds(TQueryData const & qData, LeftOverlap const) 
//{
//    return getVReadIds(qData);
//}
//
//template<typename TQueryData>
//StringSet<CharString> const & getReadIds(TQueryData const & qData, RightOverlap const) 
//{
//    return getVDJReadIds(qData);
//}

#endif
