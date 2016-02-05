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
bool belowErrRate(TSequence const & seqA, TSequence const & seqB, double const   maxErrRate)
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
    return belowErrRate(recA.fwSeq, recB.fwSeq, options.bcClustMaxErrRate) &&
        belowErrRate(recA.revSeq, recB.revSeq, options.bcClustMaxErrRate);
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

inline std::vector<FastqMultiRecord<SingleEnd>*> listFastqMultiRecords(FastqMultiRecordCollection<SingleEnd> & collection,
        String<Dna5> const & bcSeq)
{
    typedef FastqMultiRecordCollection<SingleEnd> TColl;
    typedef TColl::TBcMap TBcMap;
    typedef TColl::TSeqMap TSeqMap;
    typedef std::vector<FastqMultiRecord<SingleEnd>*> TMultiVector;

    TMultiVector result;

    TBcMap & bcMap = collection.bcMap;
    TBcMap::const_iterator bcIt = bcMap.find(bcSeq);
    if (bcIt == bcMap.end())
        return result;

    TSeqMap const & seqMap = bcIt->second;
    for (TSeqMap::const_iterator seqIt = seqMap.begin(); seqIt != seqMap.end(); ++seqIt)
    {
        FastqMultiRecord<SingleEnd> & rec = getMultiRecord(collection, seqIt->second);
        if (!rec.ids.empty())
            result.push_back(&rec);
    }
    return result;
}

inline std::vector<FastqMultiRecord<PairedEnd>*> listFastqMultiRecords(FastqMultiRecordCollection<PairedEnd> & collection,
        String<Dna5> const & bcSeq)
{
    typedef FastqMultiRecordCollection<PairedEnd> TColl;
    typedef TColl::TBcMap TBcMap;
    typedef TColl::TFwSeqMap TFwSeqMap;
    typedef TColl::TRevSeqMap TRevSeqMap;
    typedef std::vector<FastqMultiRecord<PairedEnd>*> TMultiVector;

    TMultiVector result;

    TBcMap & bcMap = collection.bcMap;
    TBcMap::const_iterator bcIt = bcMap.find(bcSeq);
    if (bcIt == bcMap.end())
        return result;

    TFwSeqMap const & fwSeqMap = bcIt->second;
    for (TFwSeqMap::const_iterator fwSeqIt = fwSeqMap.begin(); fwSeqIt != fwSeqMap.end(); ++fwSeqIt)
    {
        TRevSeqMap const & revSeqMap = fwSeqIt->second;
        for (TRevSeqMap::const_iterator revSeqIt = revSeqMap.begin(); revSeqIt != revSeqMap.end(); ++revSeqIt)
        {
            FastqMultiRecord<PairedEnd> & rec = getMultiRecord(collection, revSeqIt->second);
            if (!rec.ids.empty())
                result.push_back(&rec);
        }
    }
    return result;
}


/**
 * Collects all FastqMultiRecords per barcode ordered by size for clustering
 * @special Single end
 */
template<typename TSequencingSpec>
std::vector<std::vector<FastqMultiRecord<TSequencingSpec>*> > makeBarcodeCorrectionVector(FastqMultiRecordCollection<TSequencingSpec> & collection)
{
    typedef FastqMultiRecordCollection<TSequencingSpec> TColl;
    typedef typename TColl::TBcMap TBcMap;
    typedef std::vector<FastqMultiRecord<TSequencingSpec>*> TMultiVector;

    std::vector<TMultiVector> mVectors;
    TBcMap & bcMap = collection.bcMap;
    for (typename TBcMap::iterator bcIt = bcMap.begin(); bcIt != bcMap.end(); ++bcIt)
    {
        TMultiVector mVector = listFastqMultiRecords(collection, bcIt->first);
        std::sort(mVector.begin(), mVector.end(), FastqMultiRecordSizeComparator<TSequencingSpec>());
        mVectors.push_back(mVector);
    }
    return mVectors;
}

template<typename TSequence>
bool hammingDistAtMost(TSequence const & seqA, TSequence const & seqB, unsigned maxDist)
{
    if (length(seqA) != length(seqB))
        return false;
    size_t nErrors = 0;
    for (size_t i = 0; i<length(seqA); ++i)
    {
        if (seqA[i] != seqB[i])
            ++nErrors;
        if (nErrors > maxDist)
            return false;
    }
    return true;
}

inline void printMultiRecord(FastqMultiRecord<SingleEnd> const &) {}


inline void printMultiRecord(FastqMultiRecord<PairedEnd> const & rec)
{
    std::cerr << "===== FastqMultiRecord =====\nBC =" << rec.bcSeq << "\nFW =" << rec.fwSeq << "\nREV=" << rec.revSeq;
    for (auto x : rec.ids)
        std::cerr << "\n" << x;
    std::cerr << "\n============================\n";
}

template <typename TSequencingSpec>
void joinBarcodes(FastqMultiRecordCollection<TSequencingSpec> & collection,
        String<Dna5> const & bcSeqTarget,
        String<Dna5> const & bcSeqSource)
{
    typedef FastqMultiRecordCollection<TSequencingSpec> TColl;
    typedef typename TColl::TBcMap TBcMap;

    typename TBcMap::iterator targetIt = collection.bcMap.find(bcSeqTarget);
    SEQAN_CHECK(targetIt != collection.bcMap.end(), "Please report this error.");
    typename TBcMap::iterator sourceIt= collection.bcMap.find(bcSeqSource);
    SEQAN_CHECK(sourceIt != collection.bcMap.end(), "Please report this error.");

    // Find all records that have to be re-assigned
    std::vector<FastqMultiRecord<TSequencingSpec>*> toReassign = listFastqMultiRecords(collection, bcSeqSource);

    // Re-assign them one by one
    for (FastqMultiRecord<TSequencingSpec> * oldMultiRec : toReassign)
    {
        // Insert or update record corresponding to new barcode sequence
        FastqMultiRecord<TSequencingSpec> copy = *oldMultiRec;
        copy.bcSeq = bcSeqTarget;
        mergeRecord(collection, copy);
        // Empty the old FastqMultiRecord
        oldMultiRec->ids.clear();
    }
}

template <typename TSequencingSpec>
void clusterBarcodeSequences(FastqMultiRecordCollection<TSequencingSpec> & collection,
        CdrOptions const & options)
{
    BarcodeStats stats = getBarcodeStats(collection);

    // Order barcodes by #reads
    BarcodeStats::TBcSeqs bcSeqs = stats.bcSeqs;
    std::vector<size_t> idx;
    for (size_t i=0; i<length(bcSeqs); ++i)
        idx.push_back(i);
    std::sort(idx.begin(), idx.end(), [&stats](size_t const A, size_t const B) {return stats.nReads[A] < stats.nReads[B];});

    // Handle smallest to largest
    for (size_t const minor_idx : idx)
    {
        for (std::vector<size_t>::reverse_iterator major_idx_it = idx.rbegin(); major_idx_it != idx.rend(); ++major_idx_it)
        {
            size_t const major_idx = *major_idx_it;
            if (minor_idx == major_idx)
                break;
            if (hammingDistAtMost(bcSeqs[major_idx], bcSeqs[minor_idx], options.barcodeMaxError)) {
                joinBarcodes(collection, bcSeqs[major_idx], bcSeqs[minor_idx]);
                break;
            }
        }
    }
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
