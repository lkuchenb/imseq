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

#ifndef IMSEQ_FASTQ_MULTI_RECORD_H
#define IMSEQ_FASTQ_MULTI_RECORD_H

#include <iterator>

#include "fastq_io_types.h"
#include "fastq_multi_record_types.h"
#include "reject.h"
#include "progress_bar.h"
#include "input_information.h"
#include "sequence_data_types.h"

using namespace seqan;

/*-------------------------------------------------------------------------------
 - FastqMultiRecordCollection
 -------------------------------------------------------------------------------*/

inline std::string toString(FastqMultiRecord<SingleEnd> const & rec)
{
    std::stringstream ss;
    ss << rec.ids.size() << '\t' << rec.bcSeq << '\t' << rec.seq;
    return ss.str();
}

inline std::string toString(FastqMultiRecord<PairedEnd> const & rec)
{
    std::stringstream ss;
    ss << rec.ids.size() << '\t' << rec.bcSeq << '\t' << rec.fwSeq << '\t' << rec.revSeq;
    return ss.str();
}

/**
 * Compact a FastqMultiRecordCollection, i.e. free all FastqMultiRecords that
 * have no members (ids is empty) and remove all mappings related to those.
 * @special Single end
 */
inline void compact(FastqMultiRecordCollection<SingleEnd> & collection)
{
    typedef std::list<FastqMultiRecord<SingleEnd> >::iterator TListIt;
    std::vector<FastqMultiRecord<SingleEnd>*> delPtrs;
    std::vector<TListIt> delListIts;
    // Collect list iterators pointing to empty FastqMultiRecords
    for (TListIt listIt = collection.multiRecords.begin(); listIt != collection.multiRecords.end(); ++listIt)
        if (listIt->ids.empty())
            delListIts.push_back(listIt);
    for (TListIt listIt : delListIts)
    {
        FastqMultiRecord<SingleEnd> & rec = *listIt;
        if (collection.bcMap[listIt->bcSeq].erase(&(*listIt)) != 1) throw std::runtime_error("Please report this error");
        if (collection.seqMap[listIt->seq].erase(&(*listIt)) != 1) throw std::runtime_error("Please report this error");
        collection.multiRecords.erase(listIt);
    }
}

/**
 * @special Paired end
 */
inline void compact(FastqMultiRecordCollection<PairedEnd> & collection)
{
    typedef std::list<FastqMultiRecord<PairedEnd> >::iterator TListIt;
    std::vector<FastqMultiRecord<PairedEnd>*> delPtrs;
    std::vector<TListIt> delListIts;
    // Collect list iterators pointing to empty FastqMultiRecords
    for (TListIt listIt = collection.multiRecords.begin(); listIt != collection.multiRecords.end(); ++listIt)
        if (listIt->ids.empty())
            delListIts.push_back(listIt);
    for (TListIt listIt : delListIts)
    {
        FastqMultiRecord<PairedEnd> & rec = *listIt;
        if (collection.bcMap[listIt->bcSeq].erase(&(*listIt)) != 1) throw std::runtime_error("Please report this error");
        if (collection.fwSeqMap[listIt->fwSeq].erase(&(*listIt)) != 1) throw std::runtime_error("Please report this error");
        if (collection.revSeqMap[listIt->revSeq].erase(&(*listIt)) != 1) throw std::runtime_error("Please report this error");
        collection.multiRecords.erase(listIt);
    }
}

/**
 * Make a FastqRecord from a FastqMultiRecord, with an empty ID and default qualities
 */
inline FastqRecord<SingleEnd> toFastqRecordSkel(FastqMultiRecord<SingleEnd> const & mRec)
{
    FastqRecord<SingleEnd> rec;
    rec.seq    = mRec.seq;
    rec.bcSeq  = mRec.bcSeq;
    return rec;
}

inline FastqRecord<PairedEnd> toFastqRecordSkel(FastqMultiRecord<PairedEnd> const & mRec)
{
    FastqRecord<PairedEnd> rec;
    rec.fwSeq  = mRec.fwSeq;
    rec.revSeq = mRec.revSeq;
    rec.bcSeq  = mRec.bcSeq;
    return rec;
}

template<typename TSequencingSpec>
FastqMultiRecord<TSequencingSpec> & getMultiRecord(FastqMultiRecordCollection<TSequencingSpec> & coll, size_t const idx)
{
    SEQAN_CHECK(idx >= 0 && idx < coll.multiRecords.size(), "Please report this error.");
    return coll.multiRecords[idx];
}

template<typename TSequencingSpec>
FastqMultiRecord<TSequencingSpec> const & getMultiRecord(FastqMultiRecordCollection<TSequencingSpec> const & coll, size_t const idx)
{
    SEQAN_CHECK(idx >= 0 && idx < coll.multiRecords.size(), "Please report this error.");
    return coll.multiRecords[idx];
}

/**
 * Clear all data from a FastqMultiRecordCollection
 */
inline void clear(FastqMultiRecordCollection<SingleEnd> & coll) {
    coll.multiRecords.clear();
    coll.bcMap.clear();
    coll.seqMap.clear();
}

inline void clear(FastqMultiRecordCollection<PairedEnd> & coll) {
    coll.multiRecords.clear();
    coll.bcMap.clear();
    coll.fwSeqMap.clear();
    coll.revSeqMap.clear();
}

/**
 * Generates BarcodeStats given a FastqMultiRecordCollection
 *
 * @param coll The FastqMultiRecordCollection
 * @return The barcode usage stats
 */
template <typename TSequencingSpec>
BarcodeStats getBarcodeStats(FastqMultiRecordCollection<TSequencingSpec> const & coll) {
    typedef typename FastqMultiRecordCollection<TSequencingSpec>::TBcMap TBcMap;
    TBcMap const & bcMap = coll.bcMap;
    typedef FastqMultiRecord<TSequencingSpec> TMRec;
    BarcodeStats stats;

    for (typename TBcMap::const_iterator bcIt = bcMap.begin(); bcIt != bcMap.end(); ++bcIt) {
        uint64_t numUnique = 0;
        uint64_t numTotal = 0;
        for (TMRec const * const recPtr: bcIt->second)
        {
            FastqMultiRecord<TSequencingSpec> const & rec = *recPtr;
            uint64_t nReads = rec.ids.size();
            numTotal += nReads;
            if (nReads > 0)
                ++numUnique;
        }
        if (numUnique > 0) {
            appendValue(stats.bcSeqs, bcIt->first);
            appendValue(stats.nReads, numTotal);
            appendValue(stats.nUniqueReads, numUnique);
            stats.nTotalUniqueReads += numUnique;
            stats.nTotalReads += numTotal;
        }
    }

    return stats;
}

inline void writeBarcodeStats(BarcodeStats const & bcStats, std::string const & path)
{
    std::ofstream ofs(path);
    if (!ofs.good())
    {
        std::cerr << "\n[ERROR] Cannot open file '" << path << "' to write barcode stats\n";
        std::exit(1);
    }
    ofs << "BarcodeSeq\tnReads\tnUniqueReads\n";
    for (size_t i=0; i<length(bcStats.bcSeqs); ++i)
        ofs << bcStats.bcSeqs[i] << '\t' << bcStats.nReads[i] << '\t' << bcStats.nUniqueReads[i] << '\n';
    ofs.close();
}

/**
 * Returns the position of a FastqMultiRecord in a FastqMultiRecordCollection
 * given the sequences
 *
 * @returns Position of the corresponding FastqMultiRecord in the collection if
 *          there is a match, NO_MATCH otherwise.
 * @special Single end
 */
inline FastqMultiRecord<SingleEnd>* getMultiRecordPtr(FastqMultiRecordCollection<SingleEnd> const & collection,
        FastqRecord<SingleEnd>::TSequence const & bcSeq,
        FastqRecord<SingleEnd>::TSequence const & seq)
{
    typedef FastqMultiRecordCollection<SingleEnd> TColl;
    typedef typename TColl::TMRec *               TPtr;

    // Get the index sets for bc and seq
    TColl::TBcMap::const_iterator bcSetIt = collection.bcMap.find(bcSeq);
    TColl::TSeqMap::const_iterator seqSetIt = collection.seqMap.find(seq);
    if (bcSetIt == collection.bcMap.end() || seqSetIt == collection.seqMap.end())
        return nullptr;

    // Intersect the sets
    std::set<TPtr> intersectSet;
    std::set_intersection(bcSetIt->second.begin(), bcSetIt->second.end(), seqSetIt->second.begin(), seqSetIt->second.end(), std::inserter(intersectSet, intersectSet.end()));
    if (intersectSet.empty())
        return nullptr;
    if (intersectSet.size() > 1)
        throw std::runtime_error("For SingleEnd collections, (bcSeq,seq) tuples should be unique!");

    return *intersectSet.begin();
}

/**
 * @special Paired end
 */
inline FastqMultiRecord<PairedEnd>* getMultiRecordPtr(FastqMultiRecordCollection<PairedEnd> const & collection,
        FastqRecord<PairedEnd>::TSequence const & bcSeq,
        FastqRecord<PairedEnd>::TSequence const & fwSeq,
        FastqRecord<PairedEnd>::TSequence const & revSeq)
{
    typedef FastqMultiRecordCollection<PairedEnd> TColl;
    typedef typename TColl::TMRec *               TPtr;

    // Get the index sets for bc, fwSeq and revSeq
    TColl::TBcMap::const_iterator bcSetIt = collection.bcMap.find(bcSeq);
    TColl::TFwSeqMap::const_iterator fwSeqSetIt = collection.fwSeqMap.find(fwSeq);
    TColl::TRevSeqMap::const_iterator revSeqSetIt = collection.revSeqMap.find(revSeq);
    if (bcSetIt == collection.bcMap.end() || fwSeqSetIt == collection.fwSeqMap.end() || revSeqSetIt == collection.revSeqMap.end())
        return nullptr;

    // Intersect the sets
    std::set<TPtr> intersectSetTmp;
    std::set_intersection(bcSetIt->second.begin(), bcSetIt->second.end(), fwSeqSetIt->second.begin(), fwSeqSetIt->second.end(), std::inserter(intersectSetTmp, intersectSetTmp.end()));
    if (intersectSetTmp.empty())
        return nullptr;
    std::set<TPtr> intersectSet;
    std::set_intersection(intersectSetTmp.begin(), intersectSetTmp.end(), revSeqSetIt->second.begin(), revSeqSetIt->second.end(), std::inserter(intersectSet, intersectSet.end()));
    if (intersectSet.empty())
        return nullptr;
    if (intersectSet.size() > 1)
        throw std::runtime_error("For PairedEnd collections, (bcSeq,fwSeq,revSeq) tuples should be unique!");

    return *intersectSet.begin();
}

inline FastqMultiRecord<SingleEnd>* getMultiRecordPtr(FastqMultiRecordCollection<SingleEnd> const & collection,
        FastqRecord<SingleEnd> const & rec) {
    return getMultiRecordPtr(collection, rec.bcSeq, rec.seq);
}

inline FastqMultiRecord<PairedEnd>* getMultiRecordPtr(FastqMultiRecordCollection<PairedEnd> const & collection,
        FastqRecord<PairedEnd> const & rec) {
    return getMultiRecordPtr(collection, rec.bcSeq, rec.fwSeq, rec.revSeq);
}

template <typename TSequencingSpec>
FastqMultiRecord<TSequencingSpec> _generic_newMultiRecord(FastqRecord<TSequencingSpec> const & record) {
    FastqMultiRecord<TSequencingSpec> multiRecord;
    multiRecord.bcSeq = record.bcSeq;
    multiRecord.ids.insert(record.id);
    return(multiRecord);
}

inline void updateMeanQualityValues(String<double> & targetQualities,
        uint64_t const targetWeight,
        String<double> const & newQualities,
        uint64_t const newWeight)
{
    SEQAN_CHECK((length(targetQualities)==0 && targetWeight==0) || (length(targetQualities) == length(newQualities) && length(targetQualities) > 0), "Please report this error.");
    if (length(targetQualities) == 0)
    {
        targetQualities = newQualities;
    } else {
        for (size_t i = 0; i<length(targetQualities); ++i)
            targetQualities[i] = (static_cast<double>(targetWeight) * targetQualities[i]
                    + static_cast<double>(newWeight) * newQualities[i]) / (targetWeight + newWeight);
    }
}

/**
 * Updates mean quality values based on the current mean quality values, the
 * number of elements that contributed to them (origWeight) and the new quality
 * values
 */
template<typename TQualSequence>
void updateMeanQualityValues(String<double> & qualities,
        uint64_t const origWeight,
        TQualSequence const & seq) {
    SEQAN_CHECK((length(qualities) == 0 && origWeight == 0) || (length(qualities) == length(seq) && origWeight > 0), "Please report this error");
    if (origWeight == 0)
        resize(qualities, length(seq), 0);
    for (unsigned i=0; i<length(seq); ++i)
        qualities[i] = ((qualities[i] * origWeight) + getQualityValue(seq[i])) / (origWeight+1);
}

inline void updateMeanQualityValues(FastqMultiRecord<SingleEnd> & target,
        FastqMultiRecord<SingleEnd> const & source)
{
    updateMeanQualityValues(target.qualities, target.ids.size(), source.qualities, source.ids.size());
}

inline void updateMeanQualityValues(FastqMultiRecord<PairedEnd> & target,
        FastqMultiRecord<PairedEnd> const & source)
{
    updateMeanQualityValues(target.fwQualities, target.ids.size(), source.fwQualities, source.ids.size());
    updateMeanQualityValues(target.revQualities, target.ids.size(), source.revQualities, source.ids.size());
}

/**
 * Creates a new FastqMultiRecord based on a FastqRecord
 */
inline FastqMultiRecord<SingleEnd> newMultiRecord(FastqRecord<SingleEnd> const & record) {
    FastqMultiRecord<SingleEnd> multiRecord = _generic_newMultiRecord(record);
    multiRecord.seq = record.seq;
    updateMeanQualityValues(multiRecord.qualities, 0, record.seq);
    return(multiRecord);
}

inline FastqMultiRecord<PairedEnd> newMultiRecord(FastqRecord<PairedEnd> const & record) {
    FastqMultiRecord<PairedEnd> multiRecord = _generic_newMultiRecord(record);
    multiRecord.fwSeq = record.fwSeq;
    updateMeanQualityValues(multiRecord.fwQualities, 0, record.fwSeq);
    multiRecord.revSeq = record.revSeq;
    updateMeanQualityValues(multiRecord.revQualities, 0, record.revSeq);
    return(multiRecord);
}

/**
 * Appends a FastqMultiRecord to a FastqMultiRecordCollection and updates the
 * internal maps
 */
inline FastqMultiRecord<PairedEnd> & mapMultiRecord(FastqMultiRecordCollection<PairedEnd> & collection,
        FastqMultiRecord<PairedEnd> const & multiRecord) {
    FastqMultiRecord<PairedEnd> & newRec = * collection.multiRecords.insert(collection.multiRecords.end(), multiRecord);
    collection.bcMap[multiRecord.bcSeq].insert(&newRec);
    collection.fwSeqMap[multiRecord.fwSeq].insert(&newRec);
    collection.revSeqMap[multiRecord.revSeq].insert(&newRec);
    return newRec;
}

inline FastqMultiRecord<SingleEnd> & mapMultiRecord(FastqMultiRecordCollection<SingleEnd> & collection,
        FastqMultiRecord<SingleEnd> const & multiRecord) {
    FastqMultiRecord<SingleEnd> & newRec = * collection.multiRecords.insert(collection.multiRecords.end(), multiRecord);
    collection.bcMap[multiRecord.bcSeq].insert(&newRec);
    collection.seqMap[multiRecord.seq].insert(&newRec);
    return newRec;
}

/**
 * Adds a FastqRecord to a FastqMultiRecord by adding the ID and updating the
 * median qualities. No checking for sequence identity is performed!
 *
 * @param multiRecord The FastqMultiRecord object to modify
 * @param      record The FastqRecord to add to the FastqMultiRecord
 */
inline void updateMultiRecord(FastqMultiRecord<SingleEnd> & multiRecord,
        FastqRecord<SingleEnd> const & record) {
    uint64_t old_size = multiRecord.ids.size();
    multiRecord.ids.insert(record.id);
    SEQAN_CHECK(old_size < multiRecord.ids.size(), "Please report this error");
    updateMeanQualityValues(multiRecord.qualities, old_size, record.seq);
}

inline void updateMultiRecord(FastqMultiRecord<PairedEnd> & multiRecord,
        FastqRecord<PairedEnd> const & record) {
    uint64_t old_size = multiRecord.ids.size();
    multiRecord.ids.insert(record.id);
    SEQAN_CHECK(old_size < multiRecord.ids.size(), "Please report this error");
    updateMeanQualityValues(multiRecord.fwQualities, old_size, record.fwSeq);
    updateMeanQualityValues(multiRecord.revQualities, old_size, record.revSeq);
}

/**
 * Find the FastqMultiRecord that contains a certain FastqRecord within a
 * FastqMultiRecordCollection. If the specified FastqRecord is not yet in the
 * FastqMultiRecordCollection and insert=true is specified, add the record to
 * the collection. Matching does not consider FASTQ ID, the 'insert' feature
 * does.
 *
 * @param multiRecord The output FastqMultiRecord object to write the maching
 *                    multi-record to if one was found or inserted.
 * @param  collection The FastqMultiRecordCollection to search in / modify
 * @param      record The FastqRecord to look for / insert
 * @param      insert true = insert if no match, false = don't insert. Default: false;
 *
 * @return            A pointer to the found record, NULL if there was no
 *                    matching record.
 */
template <typename TSequencingSpec>
FastqMultiRecord<TSequencingSpec> * findContainingMultiRecord(FastqMultiRecordCollection<TSequencingSpec> & collection,
        FastqRecord<TSequencingSpec> const & record,
        bool insert = false)
{
    typedef FastqMultiRecordCollection<TSequencingSpec>      TColl;
    typedef typename TColl::TMRec                            TMRec;
    typedef typename FastqMultiRecord<TSequencingSpec>::TIds TIds;

    TMRec * recPtr = getMultiRecordPtr(collection, record);
    if (recPtr != nullptr) {
        TMRec & oldMultiRecord = *recPtr;
        TIds & ids = oldMultiRecord.ids;
        if (ids.find(record.id) != ids.end()) {
            return & oldMultiRecord;
        } else {
            if (!insert)
                return & oldMultiRecord;
            updateMultiRecord(oldMultiRecord, record);
            return & oldMultiRecord;
        }
    } else {
        if (!insert)
            return NULL;
        return &(mapMultiRecord(collection, newMultiRecord(record)));
    }
}

/**
 * Merge a FastqMultiRecord into a FastqMultiRecordCollection. If there already
 * exists a FastqMultiRecord with the very same sequence specifications, the
 * mean qualities are updated and the FASTQ_IDS are added. If no
 * FastqMultiRecord matches, the passed FastqMultiRecord is inserted and mapped
 * accordingly.
 */
template<typename TSequencingSpec>
FastqMultiRecord<TSequencingSpec> & mergeRecord(FastqMultiRecordCollection<TSequencingSpec> & collection,
        FastqMultiRecord<TSequencingSpec> const & rec)
{
    FastqMultiRecord<TSequencingSpec> * existingRecPtr = findContainingMultiRecord(collection, toFastqRecordSkel(rec), false);
    if (existingRecPtr != NULL)
    {
        FastqMultiRecord<TSequencingSpec> & existingRec = *existingRecPtr;
        updateMeanQualityValues(existingRec, rec);
        existingRec.ids.insert(rec.ids.begin(), rec.ids.end());
        existingRec.bcSeqHistory.insert(rec.bcSeqHistory.begin(), rec.bcSeqHistory.end());
        return existingRec;
    }
    return mapMultiRecord(collection, rec);
}

/**
 * Read all or at most 'counts' records from the input streams and perform
 * barcode splitting if specified.
 * @param qData The QueryData object to write to.
 * @param inStreams The SeqInputStreams to read from
 * @param options User specified options
 * @param count the maximum number of records to read. If set to '0', read until streams are exhausted
 * @return 'true' if the input streams are not yet exhausted, 'false' otherwise
 */
template <typename TSequencingSpec>
bool readRecords(FastqMultiRecordCollection<TSequencingSpec> & collection,
        InputInformation & ii,
        String<RejectEvent> & rejectEvents,
        SeqInputStreams<TSequencingSpec> & inStreams,
        CdrOptions const & options,
        unsigned count = 0)
{
    // ============================================================================
    // Initialize InputInformation
    // ============================================================================
    ii.maxReadLength  = 0;
    ii.minReadLength  = std::numeric_limits<unsigned>::max();
    ii.totalReadCount = 0;

    // ============================================================================
    // Initialize Progress bar
    // ============================================================================

    ProgressBar * progBar = NULL;
    if (inStreams.totalInBytes > 0)
        progBar = new ProgressBar(std::cerr, inStreams.totalInBytes, 100, "      ");

    // ============================================================================
    // Fill collection
    // ============================================================================

    clear(collection);
    FastqRecord<TSequencingSpec> rec;
    uint64_t blockBytes = 0;
    while (!inStreamsAtEnd(inStreams)) {
        bool tsfb = false;
        if (count > 0 && ii.totalReadCount == count) {
            progBar->clear();
            return !inStreamsAtEnd(inStreams);
        }
        if (options.barcodeLength > 0) {
            if (!readRecord(rec, inStreams, options.barcodeVDJRead, options.barcodeLength)) {
                tsfb = true;
            }
        } else {
            readRecord(rec, inStreams);
        }

        // Count the read record
        ++ii.totalReadCount;
        blockBytes += approxSizeInBytes(rec);
        if (ii.totalReadCount % 1234 == 0  && progBar != NULL) {
            progBar->updateAndPrint(blockBytes);
            blockBytes = 0;
        }

        // Truncate record if requested
        if (options.trunkReads != 0)
            truncate(rec, options.trunkReads);
        // FASTQ-Read QC
        RejectReason r = tsfb ? TOO_SHORT_FOR_BARCODE : qualityControl(rec, options);
        if (r == NONE) {
            // Statistics
            ii.maxReadLength = ii.maxReadLength > length(longerSeq(rec)) ? ii.maxReadLength : length(longerSeq(rec));
            ii.minReadLength = ii.minReadLength < length(shorterSeq(rec)) ? ii.minReadLength : length(shorterSeq(rec));
            // Sync orientation for further processing
            syncOrientation(rec, options);
            // Insert into collection
            findContainingMultiRecord(collection, rec, true);
        } else {
            appendValue(rejectEvents, RejectEvent(rec.id, r));
        }
    }
    progBar->clear();
    delete progBar;
    return false;
}

/**
 * Build a QueryDataCollection from a String of FastqMultiRecords
 *
 * @special Single end
 */
inline QueryDataCollection<SingleEnd> buildQDCollection(String<FastqMultiRecord<SingleEnd> const *> const & ptrs)
{
    QueryDataCollection<SingleEnd> qdc;
    for (FastqMultiRecord<SingleEnd> const * ptr : ptrs) {
        appendValue(qdc.queryData.seqs, ptr->seq);
        appendValue(qdc.queryData.avgQVals, ptr->qualities);
    }
    return qdc;
}

/**
 * @special Paired end
 */
inline QueryDataCollection<PairedEnd> buildQDCollection(String<FastqMultiRecord<PairedEnd> const *> const & ptrs)
{
    QueryDataCollection<PairedEnd> qdc;

    for (size_t i=0; i<length(ptrs); ++i)
    {
        FastqMultiRecord<PairedEnd> const * ptr = ptrs[i];
        if (!empty(ptr->fwSeq))
        {
            appendValue(qdc.pairedQueryData.fwSeqs, ptr->fwSeq);
            appendValue(qdc.pairedQueryData.revSeqs, ptr->revSeq);
            appendValue(qdc.pairedQueryData.fwAvgQVals, ptr->fwQualities);
            appendValue(qdc.pairedQueryData.revAvgQVals, ptr->revQualities);
        } else {
            appendValue(qdc.singleQueryData.seqs, ptr->revSeq);
            appendValue(qdc.singleQueryData.avgQVals, ptr->revQualities);
            appendValue(qdc.sePositions, i);
        }
    }
    return qdc;
}

#endif
