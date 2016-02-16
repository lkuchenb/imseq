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
// Reading FASTA / FASTQ files
// ============================================================================

#ifndef IMSEQ_FASTQ_IO_H
#define IMSEQ_FASTQ_IO_H

#include <string>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "sequence_data.h"
#include "barcode_correction.h"
#include "reject.h"
#include "collection_utils.h"
#include "qc_basics.h"


/********************************************************************************
 * STRUCTS AND CLASSES
 ********************************************************************************/

#include "fastq_io_types.h"

/********************************************************************************
 * FUNCTIONS
 ********************************************************************************/

/**
 * Computes the size of an input file
 */
inline uint64_t computeFileSize(std::string const & path) {
    SeqFileIn seqFileIn;
    open(seqFileIn, path.c_str());
    uint64_t res = 0;
    for (; !atEnd(seqFileIn.iter); ++seqFileIn.iter)
        ++res;
    return res;
}

/**
 * Returns the longer sequence of a FastqRecord. Simply returns the sequence
 * for a single end record, the longer sequence of the two for a paired end
 * record.
 *
 * @special Single end
 */
inline FastqRecord<SingleEnd>::TSequence & longerSeq(FastqRecord<SingleEnd> & rec)
{
    return rec.seq;
}

/**
 * @special Paired end
 */
inline FastqRecord<PairedEnd>::TSequence & longerSeq(FastqRecord<PairedEnd> & rec)
{
    if (length(rec.fwSeq) > length(rec.revSeq))
        return rec.fwSeq;
    return rec.revSeq;
}

/**
 * Returns the longer sequence of a FastqRecord. Simply returns the sequence
 * for a single end record, the longer sequence of the two for a paired end
 * record.
 *
 * @special Single end
 */
inline FastqRecord<SingleEnd>::TSequence & shorterSeq(FastqRecord<SingleEnd> & rec)
{
    return rec.seq;
}

/**
 * @special Paired end
 */
inline FastqRecord<PairedEnd>::TSequence & shorterSeq(FastqRecord<PairedEnd> & rec)
{
    if (length(rec.fwSeq) < length(rec.revSeq))
        return rec.fwSeq;
    return rec.revSeq;
}

/**
 * Truncate the sequences of a FastqRecord
 *
 * @special Single end
 */
inline void truncate(FastqRecord<SingleEnd> & rec, size_t len)
{
    if (length(rec.seq) > len)
        resize(rec.seq, len);
}

/**
 * @special Paired end
 */
inline void truncate(FastqRecord<PairedEnd> & rec, size_t len)
{
    if (length(rec.fwSeq) > len)
        resize(rec.fwSeq, len);
    if (length(rec.revSeq) > len)
        resize(rec.revSeq, len);
}

/**
 * Create a string representation of a FastqRecord
 * @special Paired end implementation
 */
inline std::string toString(FastqRecord<PairedEnd> const & rec) {
    std::stringstream ss;
    ss << "BARCODE\t" << rec.bcSeq << "\tFORWARD\t" << rec.fwSeq
            << "\tREVERSE\t" << rec.revSeq;
    return ss.str();
}

/**
 * Create a string representation of a FastqRecord
 * @special Single end implementation
 */
inline std::string toString(FastqRecord<SingleEnd> const & rec) {
    std::stringstream ss;
    ss << "BARCODE\t" << rec.bcSeq << "\tREAD\t" << rec.seq;
    return ss.str();
}

/**
 * Performs the first quality control performed directly after reading a FASTQ
 * record.
 *
 * @param     rec The FASTQ record to check
 * @param options The runtime options that contain the user spec thresholds
 *
 * @return A RejectReason object
 */
template <typename TSequencingSpec>
RejectReason qualityControl(FastqRecord<TSequencingSpec> const & rec,
        CdrOptions const & options)
{
    if (contains(rec.bcSeq, 'N'))
        return N_IN_BARCODE;
    if (anyQualityBelow(rec.bcSeq, options.bcQmin))
        return LOW_QUALITY_BARCODE_BASE;
    if (averageQualityBelow(rec, options.qmin))
        return AVERAGE_QUAL_FAIL;
    return NONE;
}

inline uint64_t approxSizeInBytes(FastqRecord<SingleEnd> const & rec) {
    return 2*length(rec.seq) + 2*length(rec.bcSeq) + length(rec.id) + 6;
}

inline uint64_t approxSizeInBytes(FastqRecord<PairedEnd> const & rec) {
    return 2*length(rec.fwSeq) + 2*length(rec.revSeq) + 2*length(rec.bcSeq) + 2*length(rec.id) + 12;
}

/**
 * Read a single FASTQ record from a SeqInputStreams object
 * @special Paired end implementation
 */
inline void readRecord(FastqRecord<PairedEnd> & fastqRecord, SeqInputStreams<PairedEnd> & inStreams) {
    readRecord(fastqRecord.id, fastqRecord.revSeq, inStreams.revStream);
    readRecord(fastqRecord.id, fastqRecord.fwSeq, inStreams.fwStream);
}

/**
 * Read a single FASTQ record from a SeqInputStreams object
 * @special Single end implementation
 */
inline void readRecord(FastqRecord<SingleEnd> & fastqRecord, SeqInputStreams<SingleEnd> & inStreams) {
    readRecord(fastqRecord.id, fastqRecord.seq, inStreams.stream);
}

/**
 * Read a single FASTQ record from a SeqInputStreams object and split barcode
 * @special Single end implementation
 */
inline bool readRecord(FastqRecord<SingleEnd> & fastqRecord, SeqInputStreams<SingleEnd> & inStreams, bool const barcodeVDJRead, unsigned const barcodeLength) {
    readRecord(fastqRecord.id, fastqRecord.seq, inStreams.stream);
    return splitBarcodeSeq(fastqRecord, barcodeVDJRead, barcodeLength);
}

/**
 * Read a single FASTQ record from a SeqInputStreams object and split barcode
 * @special Paired end implementation
 */
inline bool readRecord(FastqRecord<PairedEnd> & fastqRecord, SeqInputStreams<PairedEnd> & inStreams, bool const barcodeVDJRead, unsigned const barcodeLength) {
    readRecord(fastqRecord.id, fastqRecord.revSeq, inStreams.revStream);
    readRecord(fastqRecord.id, fastqRecord.fwSeq, inStreams.fwStream);
    return splitBarcodeSeq(fastqRecord, barcodeVDJRead, barcodeLength);
}

#endif
