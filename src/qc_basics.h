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

#ifndef IMSEQ_QC_BASICS_H
#define IMSEQ_QC_BASICS_H

using namespace seqan;

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "fastq_io_types.h"

/**
 * Checks if the average quality value assigned to the values of a String
 * object is below a given threshold
 */
template <typename TValue, typename TNumeric>
bool averageQualityBelow(String<TValue> const & string, TNumeric const thresh) {
    typedef typename Iterator<String<TValue> const, Rooted>::Type TIt;
    double q = 0;
    for (TIt it = begin(string); !atEnd(it); goNext(it))
        q += getQualityValue(*it);
    q = q / length(string);
    return q < thresh;
}

/**
 * Checks if the average quality value assigned to the sequence values of a
 * FastqRecord object is below a given threshold. The barcode sequence is not
 * checked. For paired end records the sequences are checked separately and of
 * the mean is below the threshold for at least one of them the result is true.
 *
 * @special Single end
 */
template<typename TNumeric>
bool averageQualityBelow(FastqRecord<SingleEnd> & rec, TNumeric const thresh, bool dropFwSeqIfLQ = false) {
    ignoreUnusedVariableWarning(dropFwSeqIfLQ);
    return averageQualityBelow(rec.seq, thresh);
}

/**
 * @special Paired end
 */
template<typename TNumeric>
bool averageQualityBelow(FastqRecord<PairedEnd> & rec, TNumeric const thresh, bool dropFwSeqIfLQ = false) {
    if (averageQualityBelow(rec.fwSeq, thresh))
    {
        if (!dropFwSeqIfLQ)
            return true;
        else if (!averageQualityBelow(rec.revSeq, thresh))
        {
            rec.fwSeq = "";
            return false;
        }
    }
    return averageQualityBelow(rec.revSeq, thresh);
}

/**
 * Checks if the quality value assigned to any of the values of a String is
 * below a given threshold
 */
template <typename TValue, typename TNumeric>
bool anyQualityBelow(String<TValue> const & string, TNumeric const thresh) {
    typedef typename Iterator<String<TValue> const, Rooted>::Type TIt;
    for (TIt it = begin(string); !atEnd(it); goNext(it))
        if (getQualityValue(*it) < thresh)
            return true;
    return false;
}

#endif
