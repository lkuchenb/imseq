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

#ifndef IMSEQ_UNIT_TESTS_IMSEQ_QC_BASICS_H
#define IMSEQ_UNIT_TESTS_IMSEQ_QC_BASICS_H

#include "../src/qc_basics.h"

SEQAN_DEFINE_TEST(unit_tests_imseq_qc_basics_averageQualityBelow_string)
{

    {
        String<Dna5Q> str = "ACGTTAAG";
        for (auto & s : str)
            assignQualityValue(s, 30);
        SEQAN_ASSERT(!averageQualityBelow(str, 30));
        
        assignQualityValue(str[2], 28);
        SEQAN_ASSERT(averageQualityBelow(str, 30));
        SEQAN_ASSERT(!averageQualityBelow(str, 29));
        str[2] = 'N';
        SEQAN_ASSERT(averageQualityBelow(str, 29));
    }
}


SEQAN_DEFINE_TEST(unit_tests_imseq_qc_basics_anyQualityBelow_string)
{

    {
        String<Dna5Q> str = "ACGTTAAG";
        for (auto & s : str)
            assignQualityValue(s, 30);
        SEQAN_ASSERT(!anyQualityBelow(str, 30));
        SEQAN_ASSERT(anyQualityBelow(str, 31));
        
        assignQualityValue(str[2], 28);
        SEQAN_ASSERT(anyQualityBelow(str, 30));
        str[2] = 'N';
        SEQAN_ASSERT(anyQualityBelow(str, 1));
    }
}

SEQAN_DEFINE_TEST(unit_tests_imseq_qc_basics_averageQualityBelow_fastqRecordSingleEnd)
{
    {
        FastqRecord<SingleEnd> rec;
        rec.seq = "ACGTTAAG";
        for (auto & s : rec.seq)
            assignQualityValue(s, 30);
        SEQAN_ASSERT(!averageQualityBelow(rec, 30));
        
        assignQualityValue(rec.seq[2], 28);
        SEQAN_ASSERT(averageQualityBelow(rec, 30));
        SEQAN_ASSERT(!averageQualityBelow(rec, 29));
        rec.seq[2] = 'N';
        SEQAN_ASSERT(averageQualityBelow(rec, 29));


    }
}

SEQAN_DEFINE_TEST(unit_tests_imseq_qc_basics_averageQualityBelow_fastqRecordPairedEnd)
{
    {
        FastqRecord<PairedEnd> rec;
        rec.fwSeq = "ACGTTAAG";
        rec.revSeq= "TGTGGCTA";
        for (auto & s : rec.fwSeq)
            assignQualityValue(s, 30);
        for (auto & s : rec.revSeq)
            assignQualityValue(s, 30);
        SEQAN_ASSERT(!averageQualityBelow(rec, 30));
        
        assignQualityValue(rec.fwSeq[2], 28);
        SEQAN_ASSERT(averageQualityBelow(rec, 30));
        SEQAN_ASSERT(!averageQualityBelow(rec, 29));
        rec.revSeq[2] = 'N';
        SEQAN_ASSERT(averageQualityBelow(rec, 29));
    }
}

#endif
