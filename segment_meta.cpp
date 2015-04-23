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

#include "segment_meta.h"

/*
 * Constructs a SegmentMeta object from a given sequence identifier
 * Returns TRUE if the sequence identifier could successfully be parsed
 * and FALSE otherwise.
 */
bool parseMetaInformation(SegmentMeta& segMeta, seqan::CharString const & idString) {
    std::stringstream s(std::string(toCString(idString)));
    std::string text;

    // (1) Parse the gene name
    if (!std::getline(s, text, '|'))
        return false;
    segMeta.geneName = text;

    // (2) Parse the segment type
    if (!std::getline(s, text, '|'))
        return false;
    segMeta.segType = text;

    // (3) Parse the segment id
    if (!std::getline(s, text, '|'))
        return false;
    segMeta.segId = text;

    // (4) Parse the allel
    if (!std::getline(s, text, '|'))
        return false;
    segMeta.allel = text;

    // (5) Parse the motif position
    if (!std::getline(s, text, '|'))
        return false;
    char *t = new char[text.size()+1];
    std::strcpy(t, text.c_str());
    segMeta.motifPos = static_cast<unsigned>(std::atoi(t));
    delete[] t;

    return true;
}

