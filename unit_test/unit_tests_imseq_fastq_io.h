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

#ifndef IMSEQ_UNIT_TESTS_IMSEQ_FASTQ_IO_H
#define IMSEQ_UNIT_TESTS_IMSEQ_FASTQ_IO_H

#include "../src/fastq_io.h"

SEQAN_DEFINE_TEST(unit_tests_imseq_fastq_io_qualityControl)
{
    {
        CdrOptions options;
        options.qmin = 30;
        options.bcQmin = 30;
        options.singleEndFallback = false;

        FastqRecord<PairedEnd> rec;
        rec.fwSeq = "CGATGACTGAGCTAGACTGCAC";
        rec.revSeq = "GCATTAGCGTTTCA";
        rec.bcSeq = "TCGAC";
        for (auto & x : rec.fwSeq)
            assignQualityValue(x,30);
        for (auto & x : rec.revSeq)
            assignQualityValue(x,30);
        for (auto & x : rec.bcSeq)
            assignQualityValue(x,30);

        SEQAN_ASSERT_EQ(qualityControl(rec, options), NONE);
        assignQualityValue(rec.fwSeq[3], 25);
        SEQAN_ASSERT_EQ(qualityControl(rec, options), AVERAGE_QUAL_FAIL);
        options.singleEndFallback = true;
        SEQAN_ASSERT_EQ(qualityControl(rec, options), NONE);
        assignQualityValue(rec.bcSeq[2], 29);
        SEQAN_ASSERT_EQ(qualityControl(rec, options), LOW_QUALITY_BARCODE_BASE);
        rec.bcSeq[4] = 'N';
        SEQAN_ASSERT_EQ(qualityControl(rec, options), N_IN_BARCODE);
    }
}

SEQAN_DEFINE_TEST(unit_tests_imseq_fastq_io_approxSizeInBytes)
{
    {
        FastqRecord<SingleEnd> seRec;
        seRec.seq = "CGTACGTACGTACGT";
        seRec.id = "M01437:163:000000000-AAAAAA:x:x:x";
        SEQAN_ASSERT_EQ(approxSizeInBytes(seRec), 69u);
        seRec.bcSeq = "ACGGG";
        SEQAN_ASSERT_EQ(approxSizeInBytes(seRec), 79u);

        FastqRecord<PairedEnd> peRec;
        peRec.fwSeq = "CGTACGTACGTACGT";
        peRec.revSeq= "GCGTATTAATTCGCAAAGTC";
        peRec.id = "M01437:163:000000000-AAAAAA:x:x:x";
        SEQAN_ASSERT_EQ(approxSizeInBytes(peRec), 148u);
        peRec.bcSeq = "ACGGG";
        SEQAN_ASSERT_EQ(approxSizeInBytes(peRec), 158u);
    }
}

#endif
