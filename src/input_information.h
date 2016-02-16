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

#ifndef IMSEQ_INPUT_INFORMATION_H
#define IMSEQ_INPUT_INFORMATION_H

struct InputInformation {
    uint64_t totalReadCount;
    unsigned maxReadLength;
    unsigned minReadLength;
    InputInformation() : totalReadCount(0), maxReadLength(0), minReadLength(0) {}
};

template<typename TCdrInput>
bool getInputInformation(InputInformation & inputInformation, TCdrInput & input)
{

    std::cerr << "===== Preprocessing" << std::endl;

    inputInformation.totalReadCount = 0;
    inputInformation.maxReadLength = 0;
    inputInformation.minReadLength = -1u;

    seqan::SeqFileIn & seqStream = getInfoStream(input);

    while (!atEnd(seqStream)) {
        CharString id;
        String<Dna5Q> seq;

        try {
            readRecord(id, seq, seqStream);
        } catch (Exception const & e) {
            std::cerr << "[ERR] Error reading FASTQ file: " << e.what() << std::endl;
            return false;
        }

        unsigned rlen = length(seq);

        if (rlen > inputInformation.maxReadLength)
            inputInformation.maxReadLength = rlen;
        if (rlen < inputInformation.minReadLength)
            inputInformation.minReadLength = rlen;

        ++inputInformation.totalReadCount;
    }

    if (inputInformation.totalReadCount==0) {
            std::cerr << "[ERR] Error reading FASTQ file!" << std::endl;
            return false;
    }

    resetStreams(input);

    return true;
}

#endif
