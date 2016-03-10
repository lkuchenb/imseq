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

#ifndef IMSEQ_SEQUENCE_DATA_TYPES_H
#define IMSEQ_SEQUENCE_DATA_TYPES_H

#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

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

typedef String<Dna5>           TQueryDataSequence;

template<typename T>
struct QueryData {
    typedef String<Dna5>	TSequence;
};

template<>
struct QueryData<SingleEnd> : QueryData<void> {
    StringSet<TQueryDataSequence> seqs;
    StringSet<String<double> >    avgQVals;
};

template<>
struct QueryData<PairedEnd> : QueryData<void> {
    StringSet<TQueryDataSequence> fwSeqs, revSeqs;
    StringSet<String<double> >    fwAvgQVals,revAvgQVals;
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
    String<size_t> sePositions;
    QueryData<PairedEnd> pairedQueryData;
    QueryData<SingleEnd> singleQueryData;
};

#endif
