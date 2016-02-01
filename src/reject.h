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
    MOTIF_QUALITY_TOO_LOW    = 3,
    NONSENSE_IN_CDR3         = 4,
    OUT_OF_READING_FRAME     = 5,
    SEGMENT_MATCH_FAILED     = 6,
    BROKEN_CDR_BOUNDARIES    = 7,
    LOW_QUALITY_BASE_IN_CDR  = 8,
    TOO_SHORT_FOR_BARCODE    = 9,
    LOW_QUALITY_BARCODE_BASE = 10,
    N_IN_BARCODE             = 11
};

const CharString _CDRREJECTS[] = {"NONE","AVERAGE_QUAL_FAIL","MOTIF_AMBIGUOUS","MOTIF_QUALITY_TOO_LOW","NONSENSE_IN_CDR3","OUT_OF_READING_FRAME","SEGMENT_MATCH_FAILED","BROKEN_CDR_BOUNDARIES","LOW_QUALITY_BASE_IN_CDR","TOO_SHORT_FOR_BARCODE","LOW_QUALITY_BARCODE_BASE","N_IN_BARCODE"};

struct RejectEvent {
    CharString id;
    RejectReason reason;
    RejectEvent(CharString _id, RejectReason _reason) : id(_id), reason(_reason) {}
};

#endif
