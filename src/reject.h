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

#ifndef IMSEQ_REJECT_H
#define IMSEQ_REJECT_H

#include <seqan/basic.h>
#include <seqan/sequence.h>

enum RejectReason {
    NONE                     = 0,
    AVERAGE_QUAL_FAIL        = 1,
    MOTIF_AMBIGUOUS          = 2,
    NONSENSE_IN_CDR3         = 3,
    OUT_OF_READING_FRAME     = 4,
    SEGMENT_MATCH_FAILED     = 5,
    BROKEN_CDR_BOUNDARIES    = 6,
    TOO_SHORT_FOR_BARCODE    = 7,
    LOW_QUALITY_BARCODE_BASE = 8,
    N_IN_BARCODE             = 9,
    READ_TOO_SHORT           = 10,
    CDR3_TOO_SHORT           = 11
};

const CharString _CDRREJECTS[] = {"NONE","AVERAGE_QUAL_FAIL","MOTIF_AMBIGUOUS","NONSENSE_IN_CDR3","OUT_OF_READING_FRAME","SEGMENT_MATCH_FAILED","BROKEN_CDR_BOUNDARIES","TOO_SHORT_FOR_BARCODE","LOW_QUALITY_BARCODE_BASE","N_IN_BARCODE","READ_TOO_SHORT","CDR3_TOO_SHORT"};

struct RejectEvent {
    CharString id;
    RejectReason reason;
    RejectEvent(CharString _id, RejectReason _reason) : id(_id), reason(_reason) {}
};

#endif
