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
#include "thread_pool.h"

using namespace seqan;

/**
 * Separate sequence and barcode which is a prefix of the sequence
 * @param       seq     The original sequence to be modified
 * @param     bcSeq     The object to store the barcode sequence in
 * @param barcodeLength Length of the barcode
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

template<typename TSequencingSpec>
struct FastqMultiRecordSizeComparator{

    bool operator() (FastqMultiRecord<TSequencingSpec> const * const lhs,
            FastqMultiRecord<TSequencingSpec> const * const rhs) const
    {
        return lhs->ids.size() > rhs->ids.size();
    }

};

/**
 * Merge identically barcoded FastqMultiRecord given in a descendingly sorted
 * vector. The resulting vector contains FastqMultiRecord with empty ID sets
 * for those records that were reassigned to another record.
 * THREADSAFE
 */
template<typename TSequencingSpec>
void mergeIdenticallyBarcoded(std::vector<FastqMultiRecord<TSequencingSpec>*> & mVector,
        CdrOptions const & options,
        ProgressBar * progBar = NULL)
{
    typedef FastqMultiRecord<TSequencingSpec> TMRec;
    typedef std::vector<FastqMultiRecord<TSequencingSpec>*> TMultiVector;

    uint64_t cnt = 0;
    for (typename TMultiVector::reverse_iterator currentIt = mVector.rbegin(); currentIt != mVector.rend(); ++currentIt)
    {
        TMRec & current = **currentIt;
        for (typename TMultiVector::iterator largerIt = mVector.begin(); largerIt != mVector.end(); ++largerIt)
        {
            TMRec & larger = **largerIt;
            // We have reached the zone of same sized FastqMultiRecords that
            // includes the current one. Abort here, we don't want to cluster
            // to a FastqMultiRecord of the same size
            if ((static_cast<double>(current.ids.size()) / static_cast<double>(larger.ids.size())) >= options.bcClustMaxFreqRate) {
                break;
            }
            // Check if the clustering specs are met and if so, reassign IDs
            // and abort
            if (withinClusteringSpecs(current, larger, options))
            {
                larger.ids.insert(current.ids.begin(), current.ids.end());
                current.ids.clear();
                break;
            }
        }
        ++cnt;
        if (progBar != NULL && cnt % 12 == 0) {
            progBar->updateAndPrint(cnt);
            cnt=0;
        }
    }
    progBar->updateAndPrint(cnt);
}

/**
 * Collects all FastqMultiRecords per barcode ordered by size for clustering
 * @special Single end
 */
inline std::vector<std::vector<FastqMultiRecord<SingleEnd>*> > makeBarcodeCorrectionVector(FastqMultiRecordCollection<SingleEnd> & collection)
{
    typedef FastqMultiRecordCollection<SingleEnd> TColl;
    typedef TColl::TBcMap TBcMap;
    typedef TColl::TSeqMap TSeqMap;
    typedef std::vector<FastqMultiRecord<SingleEnd>*> TMultiVector;

    std::vector<TMultiVector> mVectors;
    TBcMap & bcMap = collection.bcMap;
    for (TBcMap::iterator bcIt = bcMap.begin(); bcIt != bcMap.end(); ++bcIt)
    {
        TSeqMap & seqMap = bcIt->second;
        TMultiVector mVector;
        for (TSeqMap::iterator seqIt = seqMap.begin(); seqIt != seqMap.end(); ++seqIt)
        {
            mVector.push_back(&(collection.multiRecords[seqIt->second]));
        }
        std::sort(mVector.begin(), mVector.end(), FastqMultiRecordSizeComparator<SingleEnd>());
        mVectors.push_back(mVector);
    }
    return mVectors;
}

/**
 * @special Paired end
 */
inline std::vector<std::vector<FastqMultiRecord<PairedEnd>*> > makeBarcodeCorrectionVector(FastqMultiRecordCollection<PairedEnd> & collection)
{
    typedef FastqMultiRecordCollection<PairedEnd> TColl;
    typedef TColl::TBcMap TBcMap;
    typedef TColl::TFwSeqMap TFwSeqMap;
    typedef TColl::TRevSeqMap TRevSeqMap;
    typedef std::vector<FastqMultiRecord<PairedEnd>*> TMultiVector;

    std::vector<TMultiVector> mVectors;
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
                mVector.push_back(&(collection.multiRecords[revSeqIt->second]));
            }
        }
        std::sort(mVector.begin(), mVector.end(), FastqMultiRecordSizeComparator<PairedEnd>());
        mVectors.push_back(mVector);
    }
    return mVectors;
}

/**
 * Performs the barcode correction given a collection of
 * FastqMultiRecordCollection and the runtime options specified by the user
 * @special Paired end
 */
template<typename TSequencingSpec>
void barcodeCorrection(FastqMultiRecordCollection<TSequencingSpec> & collection,
        CdrOptions const & options)
{
    typedef std::vector<FastqMultiRecord<TSequencingSpec>*> TMultiVector;

    // Collect the tasks
    std::vector<TMultiVector> mVectors = makeBarcodeCorrectionVector(collection);

    // Execute them in parallel
    BarcodeStats stats = getBarcodeStats(collection);
    ProgressBar progBar(std::cerr, stats.nTotalUniqueReads , 100, "      ");
    progBar.print_progress();

    { // Scope for ThreadPool which joins threads on destruct
#ifdef __WITHCDR3THREADS__
        ThreadPool threadPool(options.jobs);
#endif
        for (TMultiVector & mVector : mVectors)
#ifdef __WITHCDR3THREADS__
            threadPool.enqueue<void>([&]()
                    {
#endif
                    mergeIdenticallyBarcoded(mVector, options, & progBar);
#ifdef __WITHCDR3THREADS__
                    }
                    );
#endif
    }

    progBar.clear();
}

#endif
