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

#endif
