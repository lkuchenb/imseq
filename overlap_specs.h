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

#ifndef CDR3FINDER_OVERLAP_SPECS_H
#define CDR3FINDER_OVERLAP_SPECS_H

#include <seqan/align.h>

struct LeftOverlap_;
typedef Tag<LeftOverlap_> LeftOverlap;

struct RightOverlap_;
typedef Tag<RightOverlap_> RightOverlap;

typedef LeftOverlap VSegment;
typedef RightOverlap JSegment;

template<typename TSpec>
struct OverlapAlign {};

template<>
struct OverlapAlign<RightOverlap> {
    typedef     AlignConfig<true,false,true,true>      Type;
};

template<>
struct OverlapAlign<LeftOverlap> {
    typedef     AlignConfig<true,true,false,true>      Type;
};
#endif
