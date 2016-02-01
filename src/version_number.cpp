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

#include "version_number.h"

#ifdef BUILD_NUMBER
const std::string IMSEQ_VERSION::STRING = QUOTE(BUILD_NUMBER) == "" || QUOTE(BUILD_NUMBER) == "1" ? IMSEQ_FALLBACK_VERSION : QUOTE(BUILD_NUMBER);
#else
const std::string IMSEQ_VERSION::STRING = IMSEQ_FALLBACK_VERSION;
#endif
