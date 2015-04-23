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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_CLONE_H_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_CLONE_H_

#include <seqan/sequence.h>
#include <seqan/file.h>
#include "cdr_utils.h"

using namespace seqan;

/**
 * Represents a clonotype, identified by its CDR sequence, its J-Gene
 * and its V-Gene
 */
template<typename T>
struct Clone {
    typedef T TAlphabet;
    std::set<unsigned>  VIds;
    String<T>           cdrSeq;
    std::set<unsigned>  JIds;
};

template<typename T>
inline bool operator==(const Clone<T>& lhs, const Clone<T>& rhs)
{
    return (lhs.VIds==rhs.VIds && lhs.cdrSeq == rhs.cdrSeq && lhs.JIds == rhs.JIds);
}

template<typename T>
inline bool operator<(const Clone<T>& __x, const Clone<T>& __y) {
    return __x.JIds < __y.JIds
        || (!(__y.JIds < __x.JIds) && __x.VIds < __y.VIds)
        || (!(__y.JIds < __x.JIds) && !(__y.VIds < __x.VIds) && __x.cdrSeq < __y.cdrSeq);
}

template<typename T>
std::string toString(Clone<T> const & clone) {
    std::stringstream ss;
    ss << "[left: " << toSeparatedStringSTL(clone.VIds) << "; seq:" << clone.cdrSeq << "; right: " << toSeparatedStringSTL(clone.JIds) << "]";
    return(ss.str());
}

#endif

