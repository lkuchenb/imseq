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
 * Compact a map mapping to multiRecordPtr-IDs. If the corresponding
 * FastqMultiRecord is empty (no IDs), delete it, replace the pointer string
 * entry with NULL and unmap it.
 */
template <typename TMap, typename TSequencingSpec>
void compact(TMap & idMap, String<FastqMultiRecord<TSequencingSpec>*> & multiRecordPtrs)
{
    typedef typename TMap::const_iterator TIter;
    String<TIter> iters;
    for (TIter iter = idMap.begin(); iter!=idMap.end(); ++iter)
    {
        if (multiRecordPtrs[iter->second]->ids.empty())
        {
            appendValue(iters, iter);
            delete multiRecordPtrs[iter->second];
            multiRecordPtrs[iter->second] = NULL;
        }
    }
    for (TIter iter : iters)
        idMap.erase(iter);
}

/**
 * Compact a map that maps to a collection. If the collection is empty, erase
 * the map entry
 */
template <typename TMap>
void compact(TMap & map)
{
    typedef typename TMap::const_iterator TIter;
    String<TIter> iters;
    for (TIter iter = map.begin(); iter != map.end(); ++iter)
        if (iter->second.empty())
            appendValue(iters, iter);
    for (TIter iter : iters)
        map.erase(iter);
}

/**
 * Compact a FastqMultiRecordCollection, i.e. free all FastqMultiRecords that
 * have no members (ids is empty) and remove all mappings related to those.
 * @special Single end
 */
inline void compact(FastqMultiRecordCollection<SingleEnd> & collection)
{
    typedef FastqMultiRecordCollection<SingleEnd> TColl;
    typedef TColl::TBcMap TBcMap;
    typedef TColl::TSeqMap TSeqMap;
    // Erase ID-empty FastqMultiRecords
    for (TBcMap::iterator bcIt = collection.bcMap.begin(); bcIt != collection.bcMap.end(); ++bcIt)
        compact(bcIt->second, collection.multiRecordPtrs);
    compact(collection.bcMap);
    // Remap to skip NULL entries
    String<FastqMultiRecord<SingleEnd>*> newPtrs;
    for (TBcMap::iterator bcIt = collection.bcMap.begin(); bcIt != collection.bcMap.end(); ++bcIt)
    {
        for (TSeqMap::iterator seqIt = bcIt->second.begin(); seqIt != bcIt->second.end(); ++seqIt)
        {
            appendValue(newPtrs, collection.multiRecordPtrs[seqIt->second]);
            seqIt->second = length(newPtrs)-1;
        }
    }
    collection.multiRecordPtrs = newPtrs;
}

/**
 * @special Paired end
 */
inline void compact(FastqMultiRecordCollection<PairedEnd> & collection)
{
    typedef FastqMultiRecordCollection<PairedEnd> TColl;
    typedef TColl::TBcMap TBcMap;
    typedef TColl::TFwSeqMap TFwSeqMap;
    typedef TColl::TRevSeqMap TRevSeqMap;
    // Erase ID-empty FastqMultiRecords
    for (TBcMap::iterator bcIt = collection.bcMap.begin(); bcIt != collection.bcMap.end(); ++bcIt)
        for (TFwSeqMap::iterator fwIt = bcIt->second.begin(); fwIt != bcIt->second.end(); ++fwIt)
            compact(fwIt->second, collection.multiRecordPtrs);
    for (TBcMap::iterator bcIt = collection.bcMap.begin(); bcIt != collection.bcMap.end(); ++bcIt)
        compact(bcIt->second);
    compact(collection.bcMap);
    // Remap to skip NULL entries
    String<FastqMultiRecord<PairedEnd>*> newPtrs;
    for (TBcMap::iterator bcIt = collection.bcMap.begin(); bcIt != collection.bcMap.end(); ++bcIt)
    {
        for (TFwSeqMap::iterator fwSeqIt = bcIt->second.begin(); fwSeqIt != bcIt->second.end(); ++fwSeqIt)
        {
            for (TRevSeqMap::iterator revSeqIt = fwSeqIt->second.begin(); revSeqIt != fwSeqIt->second.end(); ++revSeqIt)
            {
                appendValue(newPtrs, collection.multiRecordPtrs[revSeqIt->second]);
                revSeqIt->second = length(newPtrs)-1;
            }
        }
    }
    collection.multiRecordPtrs = newPtrs;
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
    SEQAN_CHECK(idx >= 0 && idx < length(coll.multiRecordPtrs), "Please report this error.");
    return *coll.multiRecordPtrs[idx];
}

template<typename TSequencingSpec>
FastqMultiRecord<TSequencingSpec> const & getMultiRecord(FastqMultiRecordCollection<TSequencingSpec> const & coll, size_t const idx)
{
    SEQAN_CHECK(idx >= 0 && idx < length(coll.multiRecordPtrs), "Please report this error.");
    return *coll.multiRecordPtrs[idx];
}

/**
 * Clear all data from a FastqMultiRecordCollection
 */
template <typename TSequencingSpec>
inline void clear(FastqMultiRecordCollection<TSequencingSpec> & coll) {
    for (auto ptr : coll.multiRecordPtrs)
        delete ptr;
    clear(coll.multiRecordPtrs);
    coll.bcMap.clear();
}

/**
 * Generates BarcodeStats given a FastqMultiRecordCollection
 *
 * @special Single end
 */
inline BarcodeStats getBarcodeStats(FastqMultiRecordCollection<SingleEnd> const & coll) {
    typedef typename FastqMultiRecordCollection<SingleEnd>::TBcMap TBcMap;
    typedef typename FastqMultiRecordCollection<SingleEnd>::TSeqMap TSeqMap;
    TBcMap const & bcMap = coll.bcMap;
    BarcodeStats stats;

    for (TBcMap::const_iterator bcIt = bcMap.begin(); bcIt != bcMap.end(); ++bcIt) {
        uint64_t cnt = 0;
        TSeqMap const & seqMap = bcIt->second;
        uint64_t seqMapSize = 0;
        for (TSeqMap::const_iterator seqIt = seqMap.begin(); seqIt != seqMap.end(); ++seqIt) {
            cnt += getMultiRecord(coll, seqIt->second).ids.size();
            if (getMultiRecord(coll, seqIt->second).ids.size() > 0)
                ++seqMapSize;
        }
        if (cnt > 0) {
            appendValue(stats.bcSeqs, bcIt->first);
            appendValue(stats.nReads, cnt);
            appendValue(stats.nUniqueReads, seqMapSize);
            stats.nTotalUniqueReads += seqMapSize;
            stats.nTotalReads += cnt;
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
 * Generates BarcodeStats given a FastqMultiRecordCollection
 *
 * @special Paired end
 */
inline BarcodeStats getBarcodeStats(FastqMultiRecordCollection<PairedEnd> const & coll) {
    typedef typename FastqMultiRecordCollection<PairedEnd>::TBcMap TBcMap;
    typedef typename FastqMultiRecordCollection<PairedEnd>::TFwSeqMap TFwSeqMap;
    typedef typename FastqMultiRecordCollection<PairedEnd>::TRevSeqMap TRevSeqMap;
    TBcMap const & bcMap = coll.bcMap;
    BarcodeStats stats;

    for (TBcMap::const_iterator bcIt = bcMap.begin(); bcIt != bcMap.end(); ++bcIt) {
        uint64_t cnt = 0;
        uint64_t uCnt = 0;
        TFwSeqMap const & fwSeqMap = bcIt->second;
        for (TFwSeqMap::const_iterator fwSeqIt = fwSeqMap.begin(); fwSeqIt != fwSeqMap.end(); ++fwSeqIt) {
            TRevSeqMap const & revSeqMap = fwSeqIt->second;
            for (TRevSeqMap::const_iterator revSeqIt = revSeqMap.begin(); revSeqIt!=revSeqMap.end(); ++revSeqIt) {
                uint64_t s = getMultiRecord(coll, revSeqIt->second).ids.size();
                cnt += s;
                if (s>0)
                    ++uCnt;
            }
        }
        if (cnt > 0) {
            appendValue(stats.bcSeqs, bcIt->first);
            appendValue(stats.nReads, cnt);
            appendValue(stats.nUniqueReads, uCnt);
            stats.nTotalUniqueReads += uCnt;
            stats.nTotalReads += cnt;
        }
    }
    return stats;
}

/**
 * Returns the position of a FastqMultiRecord in a FastqMultiRecordCollection
 * given the sequences
 *
 * @returns Position of the corresponding FastqMultiRecord in the collection if
 *          there is a match, NO_MATCH otherwise.
 * @special Single end
 */
inline uint64_t findMultiRecordPosition(FastqMultiRecordCollection<SingleEnd> const & collection,
        FastqRecord<SingleEnd>::TSequence const & bcSeq,
        FastqRecord<SingleEnd>::TSequence const & seq)
{
    typedef FastqMultiRecordCollection<SingleEnd> TColl;
    TColl::TBcMap::const_iterator bcSearch = collection.bcMap.find(bcSeq);
    if (bcSearch == collection.bcMap.end())
        return TColl::NO_MATCH;
    TColl::TSeqMap::const_iterator seqSearch = bcSearch->second.find(seq);
    if (seqSearch == bcSearch->second.end())
        return TColl::NO_MATCH;
    return seqSearch->second;
}

/**
 * @special Paired end
 */
inline uint64_t findMultiRecordPosition(FastqMultiRecordCollection<PairedEnd> const & collection,
        FastqRecord<PairedEnd>::TSequence const & bcSeq,
        FastqRecord<PairedEnd>::TSequence const & fwSeq,
        FastqRecord<PairedEnd>::TSequence const & revSeq)
{
    typedef FastqMultiRecordCollection<PairedEnd> TColl;
    TColl::TBcMap::const_iterator bcSearch = collection.bcMap.find(bcSeq);
    if (bcSearch == collection.bcMap.end())
        return TColl::NO_MATCH;
    TColl::TFwSeqMap::const_iterator fwSeqSearch = bcSearch->second.find(fwSeq);
    if (fwSeqSearch == bcSearch->second.end())
        return TColl::NO_MATCH;
    TColl::TRevSeqMap::const_iterator revSeqSearch = fwSeqSearch->second.find(revSeq);
    if (revSeqSearch == fwSeqSearch->second.end())
        return TColl::NO_MATCH;
    return revSeqSearch->second;
}

inline uint64_t findMultiRecordPosition(FastqMultiRecordCollection<SingleEnd> const & collection,
        FastqRecord<SingleEnd> const & rec) {
    return findMultiRecordPosition(collection, rec.bcSeq, rec.seq);
}

inline uint64_t findMultiRecordPosition(FastqMultiRecordCollection<PairedEnd> const & collection,
        FastqRecord<PairedEnd> const & rec) {
    return findMultiRecordPosition(collection, rec.bcSeq, rec.fwSeq, rec.revSeq);
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
    appendValue(collection.multiRecordPtrs, new FastqMultiRecord<PairedEnd>(multiRecord));
    SEQAN_CHECK(collection.bcMap[multiRecord.bcSeq][multiRecord.fwSeq].find(multiRecord.revSeq)
            == collection.bcMap[multiRecord.bcSeq][multiRecord.fwSeq].end(), "Please report this error");
    size_t new_idx = length(collection.multiRecordPtrs) - 1;
    collection.bcMap[multiRecord.bcSeq][multiRecord.fwSeq][multiRecord.revSeq] = new_idx;
    return getMultiRecord(collection, new_idx);
}

inline FastqMultiRecord<SingleEnd> & mapMultiRecord(FastqMultiRecordCollection<SingleEnd> & collection,
        FastqMultiRecord<SingleEnd> const & multiRecord) {
    appendValue(collection.multiRecordPtrs, new FastqMultiRecord<SingleEnd>(multiRecord));
    SEQAN_CHECK(collection.bcMap[multiRecord.bcSeq].find(multiRecord.seq)
            == collection.bcMap[multiRecord.bcSeq].end(), "Please report this error");
    size_t new_idx = length(collection.multiRecordPtrs) - 1;
    collection.bcMap[multiRecord.bcSeq][multiRecord.seq] = new_idx;
    return getMultiRecord(collection, new_idx);
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
    typedef FastqMultiRecordCollection<TSequencingSpec> TColl;
    typedef typename FastqMultiRecord<TSequencingSpec>::TIds TIds;

    uint64_t mRecId = findMultiRecordPosition(collection, record);
    if (mRecId != TColl::NO_MATCH) {
        FastqMultiRecord<TSequencingSpec> & oldMultiRecord = getMultiRecord(collection, mRecId);
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

inline void printCollection(FastqMultiRecordCollection<SingleEnd> const & coll)
{
    for (auto bcMapElem : coll.bcMap)
    {
        std::cerr << bcMapElem.first << "\n";
        for (auto seqElem : bcMapElem.second)
        {
            std::cerr << "    " << seqElem.first << std::endl;
            for (auto recId : coll.multiRecordPtrs[seqElem.second]->ids)
            {
                std::cerr << "            " << recId << '\n';
            }
        }
    }
}

inline void printCollection(FastqMultiRecordCollection<PairedEnd> const & coll)
{
    for (auto bcMapElem : coll.bcMap)
    {
        std::cerr << bcMapElem.first << "\n";
        for (auto fwSeqElem : bcMapElem.second)
        {
            std::cerr << "    " << fwSeqElem.first << std::endl;
            for (auto revSeqElem : fwSeqElem.second)
            {
                std::cerr << "        " << revSeqElem.first << std::endl;
                for (auto recId : coll.multiRecordPtrs[revSeqElem.second]->ids)
                {
                    std::cerr << "            " << recId << '\n';
                }
            }
        }
    }
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
inline QueryDataCollection<SingleEnd> buildQDCollection(String<FastqMultiRecord<SingleEnd>*> const & ptrs)
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
inline QueryDataCollection<PairedEnd> buildQDCollection(String<FastqMultiRecord<PairedEnd>*> const & ptrs)
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
