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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_SEGMENT_META_H_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_SEGMENT_META_H_

#include <string>
#include <sstream>
#include <set>
#include <seqan/sequence.h>

using namespace seqan;

struct SegmentMeta {
    CharString  geneName, segType, segId, allel;
    unsigned motifPos;
};

bool parseMetaInformation(SegmentMeta &, CharString const &);

//! Appends a list of segment identifiers to a character string
inline void _appendGeneList(CharString& geneList, std::set<unsigned> const & idList, String<SegmentMeta> const & meta, bool mergeAllels) {
    std::set<CharString> segNames;
    for (std::set<unsigned>::const_iterator id = idList.begin(); id!=idList.end(); ++id) {
        SegmentMeta sm = meta[*id];
        CharString segName;
        append(segName, sm.segType);
        append(segName, sm.segId);
        if (mergeAllels && segNames.find(segName)!=segNames.end())
            continue;
        if (mergeAllels)
            segNames.insert(segName);
        if (!mergeAllels) {
            appendValue(segName, '*');
            append(segName, sm.allel);
        }
        append(geneList, segName);
    }
}

inline std::string getDescriptor(SegmentMeta const & sm) {
    std::string s = "[";
    s = s + toCString(sm.geneName) + toCString(sm.segType) + toCString(sm.segId) + "*" + toCString(sm.allel);
    s = s + "]";
    return s;
}

#endif
