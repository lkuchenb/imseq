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
// FILE DESCRIPTION
// ============================================================================
// Main executable for IMSEQ unit tests
// ============================================================================

#include "unit_tests_imseq_barcode_correction.h"
#include "unit_tests_imseq_fastq_io.h"
#include "unit_tests_imseq_qc_basics.h"

SEQAN_BEGIN_TESTSUITE(unit_tests_imseq)
{
    // unit_tests_imseq_barcode_correction.h
    SEQAN_CALL_TEST(unit_tests_imseq_barcode_correction_splitBarcodeSeq);
    SEQAN_CALL_TEST(unit_tests_imseq_barcode_correction_splitBarcodeSeq__FastqRecord);

    // unit_tests_imseq_fastq_io.h
    SEQAN_CALL_TEST(unit_tests_imseq_fastq_io_qualityControl);
    SEQAN_CALL_TEST(unit_tests_imseq_fastq_io_approxSizeInBytes);

    // unit_tests_imseq_qc_basics.h
    SEQAN_CALL_TEST(unit_tests_imseq_qc_basics_averageQualityBelow_string);
    SEQAN_CALL_TEST(unit_tests_imseq_qc_basics_anyQualityBelow_string);
    SEQAN_CALL_TEST(unit_tests_imseq_qc_basics_averageQualityBelow_fastqRecordSingleEnd);
    SEQAN_CALL_TEST(unit_tests_imseq_qc_basics_averageQualityBelow_fastqRecordPairedEnd);
}

SEQAN_END_TESTSUITE
