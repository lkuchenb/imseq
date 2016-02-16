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

#ifndef IMSEQ_UNIT_TESTS_IMSEQ_FASTQ_MULTI_RECORD_H
#define IMSEQ_UNIT_TESTS_IMSEQ_FASTQ_MULTI_RECORD_H

#include "../src/fastq_multi_record.h"

SEQAN_DEFINE_TEST(unit_tests_imseq_fastq_multi_record_findContainingMultiRecord_PairedEnd)
{
    {
        FastqMultiRecordCollection<PairedEnd> collection;
        SEQAN_ASSERT(empty(collection.multiRecordPtrs));
        SEQAN_ASSERT(collection.bcMap.empty());
        SEQAN_ASSERT(NULL == findContainingMultiRecord(collection, FastqRecord<PairedEnd>("READ_1", "ACTGGCATACG", "GGGGCAAAGCA", "CCAT"), false));
        FastqMultiRecord<PairedEnd> * ptr = findContainingMultiRecord(collection, FastqRecord<PairedEnd>("READ_2", "ACTGGCATACG", "GGGGCAAAGCA", "CCAT"), true);
        SEQAN_ASSERT(NULL != ptr);
        SEQAN_ASSERT(ptr == findContainingMultiRecord(collection, FastqRecord<PairedEnd>("READ_3", "ACTGGCATACG", "GGGGCAAAGCA", "CCAT"), true));
        SEQAN_ASSERT(ptr == findContainingMultiRecord(collection, FastqRecord<PairedEnd>("READ_4", "ACTGGCATACG", "GGGGCAAAGCA", "CCAT"), true));
        SEQAN_ASSERT(ptr != findContainingMultiRecord(collection, FastqRecord<PairedEnd>("READ_5", "ACTGGCATACG", "GGGGCAAAGCA", "CGAT"), true));
        SEQAN_ASSERT(ptr != findContainingMultiRecord(collection, FastqRecord<PairedEnd>("READ_6", "ACTGGCATACG", "GGGGCTAAGCA", "CGAT"), true));
        SEQAN_ASSERT(ptr != findContainingMultiRecord(collection, FastqRecord<PairedEnd>("READ_7", "ACTGGGATACG", "GGGGCTAAGCA", "CGAT"), true));
        SEQAN_ASSERT_EQ(ptr->ids.size(), 3u);
        SEQAN_ASSERT(ptr->ids.find("READ_2") != ptr->ids.end());
        SEQAN_ASSERT(ptr->ids.find("READ_3") != ptr->ids.end());
        SEQAN_ASSERT(ptr->ids.find("READ_4") != ptr->ids.end());
        SEQAN_ASSERT_EQ(collection.bcMap["CCAT"].size(), 1u);
        SEQAN_ASSERT_EQ(collection.bcMap["CGAT"].size(), 2u);
        SEQAN_ASSERT_EQ(collection.bcMap["CCAT"]["ACTGGCATACG"].size(), 1u);
        SEQAN_ASSERT_EQ(collection.bcMap["CGAT"]["ACTGGCATACG"].size(), 2u);
        SEQAN_ASSERT_EQ(collection.bcMap["CGAT"]["ACTGGGATACG"].size(), 1u);
    }
}

SEQAN_DEFINE_TEST(unit_tests_imseq_fastq_multi_record_findContainingMultiRecord_SingleEnd)
{
    {
        FastqMultiRecordCollection<SingleEnd> collection;
        SEQAN_ASSERT(empty(collection.multiRecordPtrs));
        SEQAN_ASSERT(collection.bcMap.empty());
        SEQAN_ASSERT(NULL == findContainingMultiRecord(collection, FastqRecord<SingleEnd>("READ_1", "ACTGGCATACG", "CCAT"), false));
        FastqMultiRecord<SingleEnd> * ptr = findContainingMultiRecord(collection, FastqRecord<SingleEnd>("READ_2", "ACTGGCATACG", "CCAT"), true);
        SEQAN_ASSERT(NULL != ptr);
        SEQAN_ASSERT(ptr == findContainingMultiRecord(collection, FastqRecord<SingleEnd>("READ_3", "ACTGGCATACG", "CCAT"), true));
        SEQAN_ASSERT(ptr == findContainingMultiRecord(collection, FastqRecord<SingleEnd>("READ_4", "ACTGGCATACG", "CCAT"), true));
        SEQAN_ASSERT(ptr != findContainingMultiRecord(collection, FastqRecord<SingleEnd>("READ_5", "ACTGGCATACG", "CGAT"), true));
        SEQAN_ASSERT(ptr != findContainingMultiRecord(collection, FastqRecord<SingleEnd>("READ_6", "ACTGGCATACG", "CGAT"), true));
        SEQAN_ASSERT(ptr != findContainingMultiRecord(collection, FastqRecord<SingleEnd>("READ_7", "ACTGGGATACG", "CGAT"), true));
        SEQAN_ASSERT_EQ(ptr->ids.size(), 3u);
        SEQAN_ASSERT(ptr->ids.find("READ_2") != ptr->ids.end());
        SEQAN_ASSERT(ptr->ids.find("READ_3") != ptr->ids.end());
        SEQAN_ASSERT(ptr->ids.find("READ_4") != ptr->ids.end());
        SEQAN_ASSERT_EQ(collection.bcMap["CCAT"].size(), 1u);
        SEQAN_ASSERT_EQ(collection.bcMap["CGAT"].size(), 2u);
    }
}


#endif
