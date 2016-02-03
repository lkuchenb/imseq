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

#include "fastq_multi_record_types.h"

using namespace seqan;

/*-------------------------------------------------------------------------------
 - FastqMultiRecord
 -------------------------------------------------------------------------------*/

std::hash<std::string> const DnaStringHash::hash_fun = std::hash<std::string>();

size_t DnaStringHash::operator()(String<Dna5> const & str) const
{
    size_t x = hash_fun(std::string(toCString(String<char>(str))));
    return x;
}

