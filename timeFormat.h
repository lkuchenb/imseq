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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_TIMEFORMAT_H_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_TIMEFORMAT_H_

#include<string>
#include<sstream>

// ============================================================================
// Functions
// ============================================================================

/**
 * Returns a string of the format <days>d <hours>:<minutes>:<seconds> given
 * a value in seconds
 */
std::string formatSeconds(double totalSeconds) {
    std::stringstream ss;

    unsigned days = static_cast<unsigned>(std::floor(totalSeconds/86400.0));
    unsigned hours = static_cast<unsigned>(std::floor(totalSeconds/3600.0));
    unsigned minutes = static_cast<unsigned>(std::floor(std::fmod(totalSeconds,3600.0)/60.0));
    unsigned seconds = static_cast<unsigned>(std::fmod(totalSeconds,60.0));

    ss << days << "d " << std::setfill('0') << std::setw(2) << hours << ":" 
        << std::setfill('0') << std::setw(2) << minutes << ":" 
        << std::setfill('0') << std::setw(2) << seconds;

    return(ss.str());
}

#endif  // #ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_TIMEFORMAT_H_
