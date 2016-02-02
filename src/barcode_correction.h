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
// Functions for random barcode based aggregation of read sequences
// ============================================================================

#ifndef SANDBOX_LKUCHENB_APPS_IMSEQ_BARCODE_CORRECTION_H
#define SANDBOX_LKUCHENB_APPS_IMSEQ_BARCODE_CORRECTION_H

#include <string>

#include <seqan/sequence.h>

#include "runtime_options.h"
#include "fastq_io_types.h"
#include "fastq_multi_record_types.h"
#include "fastq_multi_record.h"

using namespace seqan;

/**
 * Separate barcode and sequence
 * @param       seq The original sequence to be modified
 * @param     bcSeq The object to store the barcode sequence in
 * @param      opts Runtime options
 */
template<typename TSequence>
bool splitBarcodeSeq(TSequence & seq, TSequence & bcSeq, unsigned const barcodeLength) {
    typedef typename Infix<TSequence>::Type TInfix;

    if (barcodeLength == 0)
        return true;

    if (length(seq) < barcodeLength) {
        resize(bcSeq, 0);
        return false;
    }

    TInfix bc = infix(seq, 0, barcodeLength);
    bcSeq = bc;

    TInfix seqInfix = infix(seq, barcodeLength, length(seq));
    TSequence _seq = seqInfix;

    seq = _seq;

    return true;
}

inline bool splitBarcodeSeq(FastqRecord<SingleEnd> & rec, bool const barcodeVDJRead, unsigned const barcodeLength) {
    ignoreUnusedVariableWarning(barcodeVDJRead);
    return splitBarcodeSeq(rec.seq, rec.bcSeq, barcodeLength);
}

inline bool splitBarcodeSeq(FastqRecord<PairedEnd> & rec, bool const barcodeVDJRead, unsigned const barcodeLength) {
    if (barcodeVDJRead)
        return splitBarcodeSeq(rec.revSeq, rec.bcSeq, barcodeLength);
    else
        return splitBarcodeSeq(rec.fwSeq, rec.bcSeq, barcodeLength);
}

/**
 * Compute the error rate of two sequences using a global Myers Bit Vector
 * algorithm implementation to compute the Levenshtein distance and derive the
 * error rate from that
 */
template<typename TSequence>
bool belowErrRate(TSequence const & seqA, TSequence const & seqB, double maxErrRate)
{
    int seqLength  = std::max(length(seqA), length(seqB));
    int maxErrors  = static_cast<int>(std::floor(maxErrRate * static_cast<double>(seqLength)));
    int lenDiff    = length(seqA) > length(seqB) ? length(seqA) - length(seqB) : length(seqB) - length(seqA);
    if (lenDiff > maxErrors)
        return false;
    int score      = globalAlignmentScore(seqA, seqB, MyersBitVector());
    double errRate = static_cast<double>(-score) / static_cast<double>(seqLength);
    return errRate <= maxErrRate;
}

/**
 * Check if two FastqMultiRecord are sequence wise similar enough to be
 * clustered in the barcode clustering step according to the user specified
 * maximum error rate. Does not check for barcode identity.
 * @special Single end implementation
 */
inline bool withinClusteringSpecs(FastqMultiRecord<SingleEnd> const & recA,
        FastqMultiRecord<SingleEnd> const & recB,
        CdrOptions const & options)
{
    return belowErrRate(recA.seq, recB.seq, options.bcClustMaxErrRate);
}

/**
 * @special Paired end implementation
 */
inline bool withinClusteringSpecs(FastqMultiRecord<PairedEnd> const & recA,
        FastqMultiRecord<PairedEnd> const & recB,
        CdrOptions const & options)
{
    return belowErrRate(recA.fwSeq, recB.fwSeq, options.bcClustMaxErrRate)
        && belowErrRate(recA.revSeq, recB.revSeq, options.bcClustMaxErrRate);
}

template<typename TSequencingSpec, typename TIterator>
struct FastqMultiRecordSizeComparator{

    FastqMultiRecordCollection<TSequencingSpec> const & collection;

    bool operator() (TIterator const lhs,
            TIterator const rhs)
    {
        return collection.multiRecords[lhs->second].ids.size() > collection.multiRecords[rhs->second].ids.size();
    }

    FastqMultiRecordSizeComparator(FastqMultiRecordCollection<TSequencingSpec> const & collection_) : collection(collection_) {}

};

/**
 * Merge identically barcoded FastqMultiRecord given in a descendingly sorted
 * vector. The resulting vector contains FastqMultiRecord with empty ID sets
 * for those records that were reassigned to another record.
 */
template<typename TSequencingSpec, typename TIterator>
void mergeIdenticallyBarcoded(std::vector<TIterator> & mVector,
        FastqMultiRecordCollection<TSequencingSpec> & coll,
        CdrOptions const & options)
{
    typedef FastqMultiRecord<TSequencingSpec> TMRec;
    typedef std::vector<TIterator> TMultiVector;

    for (typename TMultiVector::reverse_iterator currentIt = mVector.rbegin(); currentIt != mVector.rend(); ++currentIt)
    {
        typename TMultiVector::reverse_iterator nextLargerIt = currentIt;
        TMRec & current = coll.multiRecords[(*currentIt)->second];
        for (++nextLargerIt; nextLargerIt != mVector.rend(); ++nextLargerIt)
        {
            TMRec & nextLarger = coll.multiRecords[(*nextLargerIt)->second];
            if (withinClusteringSpecs(current, nextLarger, options))
            {
                nextLarger.ids.insert(current.ids.begin(), current.ids.end());
                current.ids.clear();
            }
        }
    }
}

/**
 * Performs the barcode correction given a collection of
 * FastqMultiRecordCollection and the runtime options specified by the user
 * @special Paired end
 */
inline void barcodeCorrection(FastqMultiRecordCollection<PairedEnd> & collection, CdrOptions const & options)
{
    typedef FastqMultiRecordCollection<PairedEnd> TColl;
    typedef TColl::TBcMap TBcMap;
    typedef TColl::TFwSeqMap TFwSeqMap;
    typedef TColl::TRevSeqMap TRevSeqMap;
    typedef std::vector<TRevSeqMap::iterator> TMultiVector;

    TBcMap & bcMap = collection.bcMap;
    for (TBcMap::iterator bcIt = bcMap.begin(); bcIt != bcMap.end(); ++bcIt)
    {
        TFwSeqMap & fwSeqMap = bcIt->second;
        TMultiVector mVector;
        for (TFwSeqMap::iterator fwSeqIt = fwSeqMap.begin(); fwSeqIt != fwSeqMap.end(); ++fwSeqIt)
        {
            TRevSeqMap & revSeqMap = fwSeqIt->second;
            for (TRevSeqMap::iterator revSeqIt = revSeqMap.begin(); revSeqIt != revSeqMap.end(); ++revSeqIt)
            {
                mVector.push_back(revSeqIt);
            }
        }
        std::sort(mVector.begin(), mVector.end(), FastqMultiRecordSizeComparator<PairedEnd, TRevSeqMap::iterator>(collection));
        mergeIdenticallyBarcoded(mVector, collection, options);
    }
}

/**
 * @special Single end
 */
inline void barcodeCorrection(FastqMultiRecordCollection<SingleEnd> & collection, CdrOptions const & options)
{
    typedef FastqMultiRecordCollection<SingleEnd> TColl;
    typedef TColl::TBcMap TBcMap;
    typedef TColl::TSeqMap TSeqMap;
    typedef std::vector<TSeqMap::iterator> TMultiVector;

    TBcMap & bcMap = collection.bcMap;
    for (TBcMap::iterator bcIt = bcMap.begin(); bcIt != bcMap.end(); ++bcIt)
    {
        TSeqMap & seqMap = bcIt->second;
        TMultiVector mVector;
        for (TSeqMap::iterator seqIt = seqMap.begin(); seqIt != seqMap.end(); ++seqIt)
        {
            mVector.push_back(seqIt);
        }
        std::sort(mVector.begin(), mVector.end(), FastqMultiRecordSizeComparator<SingleEnd, TSeqMap::iterator>(collection));
        mergeIdenticallyBarcoded(mVector, collection, options);
    }
}

#endif
