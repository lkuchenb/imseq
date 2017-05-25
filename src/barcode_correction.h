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

/**
 * Compute the error rate of two sequences using a global Myers Bit Vector
 * algorithm implementation to compute the Levenshtein distance and derive the
 * error rate from that
 */
template<typename TSequence>
bool belowErrRate(TSequence const & seqA, TSequence const & seqB, double const  maxErrRate)
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
        double const maxErrRate)
{
    return belowErrRate(recA.seq, recB.seq, maxErrRate);
}

/**
 * @special Paired end implementation
 */
inline bool withinClusteringSpecs(FastqMultiRecord<PairedEnd> const & recA,
        FastqMultiRecord<PairedEnd> const & recB,
        double const maxErrRate)
{
    return belowErrRate(recA.fwSeq, recB.fwSeq, maxErrRate) &&
        belowErrRate(recA.revSeq, recB.revSeq, maxErrRate);
}

/// Mutex to be used in compareOneRead() when writing results
static std::mutex compareOneRead_result_write;

/**
  * Compare one read (pair) against more abundant ones for BC clustering.
  *
  * Compares one read (pair) (reference record) to every more abundant read
  * (pair) (target record) from the most frequent to the least frequent, until
  * it finds a target record to join the reference record with. If a reference
  * record is found, the method stores this in the passed clusterTodo vector
  * and returns. Otherwise the clusterTodo object remains unmodified.
  *
  * @param clusterTodo The to-cluster pairs vector to store a result in.
  * @param recs The underlying sequence records.
  * @param sortedIndices The recs indices sorted from the least to the most abundant record.
  * @param refIndex The sortedIndices index that contains the index of the reference record.
  * @param bcDelta The maximum hamming distance of two barcodes to trigger a join
  * @param seqErrRate The maximum pairwise sequence error rate for reads to be joined
  * @param maxRatio The maximum ratio refRecord / tarRecord to allow a join
  */
template<typename TSequencingSpec>
void compareOneRead(
        std::vector<std::pair<size_t,size_t> > & clusterTodo,
        std::vector<FastqMultiRecord<TSequencingSpec>*> const & sortedRecPtrs,
        size_t const refIndex,
        short const bcDelta,
        double const seqErrRate,
        double const maxRatio,
        ProgressBar & progBar)
{
    typedef FastqMultiRecord<TSequencingSpec> TMRec;
    std::vector<std::pair<size_t,size_t> > res;

    TMRec const & refRec = * sortedRecPtrs[refIndex];
    uint64_t cnt = 0;
    for (size_t tarIndex = length(sortedRecPtrs)-1; tarIndex > refIndex; --tarIndex)
    {
        TMRec const & tarRec = * sortedRecPtrs[tarIndex];
        double ratio = static_cast<double>(refRec.ids.size()) / tarRec.ids.size();
        uint64_t skipped = 0;
        // If the ratio is too high, we can stop looking at additional reads,
        // since we are iterating them in descending order
        if (ratio > maxRatio)
            skipped = tarIndex - refIndex;
        // If barcodes are too different, continue
        else if (!hammingDistAtMost(refRec.bcSeq, tarRec.bcSeq, bcDelta))
            ++cnt;
        // Check if the read sequence(s) are similar within specs
        else if (withinClusteringSpecs(refRec, tarRec, seqErrRate))
        {
            res.push_back(std::pair<size_t,size_t>(refIndex, tarIndex));
            skipped = tarIndex - refIndex;
        }

        cnt += skipped;
        if (cnt % 1234 == 0) {
            progBar.updateAndPrint(cnt);
            cnt=0;
        }
        if (skipped > 0)
            break;
    }
    progBar.updateAndPrint(cnt);
    // Store result
    std::lock_guard<std::mutex> lock(compareOneRead_result_write);
    clusterTodo.insert(clusterTodo.end(), res.begin(), res.end());
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
    typedef FastqMultiRecordCollection<TSequencingSpec> TColl;
    typedef typename TColl::TMRec                       TMRec;
    typedef std::pair<size_t,size_t>                    TPair;

    std::cerr << "  |   Sorting unique reads by frequency" << std::endl;

    // Sort the Record indices by record abundance
    std::vector<TMRec*> sortedRecPtrs;
    sortedRecPtrs.reserve(length(collection.multiRecords));
    for (TMRec & rec : collection.multiRecords)
        sortedRecPtrs.push_back(&rec);

    std::sort(sortedRecPtrs.begin(), sortedRecPtrs.end(),
            [& collection](TMRec const * a, TMRec const * b) -> bool {
            return a->ids.size() < b->ids.size();
            });

    // From least to most abundant
    std::vector<TPair> todoPairs;

    std::cerr << "  |   Pairwise UMI and read comparison" << std::endl;
    ProgressBar progBar(std::cerr, (sortedRecPtrs.size()-1)*(sortedRecPtrs.size()-1)/2, 100, "      ");
    progBar.print_progress();

    { // Scope for ThreadPool which joins threads on destruct
#ifdef __WITHCDR3THREADS__
        ThreadPool threadPool(options.jobs);
#endif
        for (size_t i = 0; i < sortedRecPtrs.size(); ++i)
        {
#ifdef __WITHCDR3THREADS__
            threadPool.enqueue<void>([i,&todoPairs,&collection,&sortedRecPtrs, &options, &progBar]()
                    {
#endif
                    compareOneRead(todoPairs, sortedRecPtrs, i, options.barcodeMaxError, options.bcClustMaxErrRate, options.bcClustMaxFreqRate, progBar);
#ifdef __WITHCDR3THREADS__
                    });
#endif
        }
    } // End ThreadPool

    progBar.clear();

    std::cerr << "  |   Merging identified pairs of reads" << std::endl;

    std::sort(todoPairs.begin(), todoPairs.end(), [](TPair const & a, TPair const & b) -> bool { return a.second < b.second; });

    for (auto & todoPair : todoPairs)
    {
        TMRec & refRec = * sortedRecPtrs[todoPair.first];
        TMRec & tarRec = * sortedRecPtrs[todoPair.second];
        tarRec.ids.insert(refRec.ids.begin(), refRec.ids.end());
        refRec.ids.clear();
        // Should be empty
        tarRec.bcSeqHistory.insert(refRec.bcSeqHistory.begin(), refRec.bcSeqHistory.end());
    }
}

#endif
