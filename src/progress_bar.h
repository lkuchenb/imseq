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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_PROGRESS_BAR_H_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_PROGRESS_BAR_H_

#include "thread_check.h"
#include "fixed_size_types.h"

#ifdef __WITHCDR3THREADS__
#include <mutex>
#endif
#include <string>
#include <ostream>
#include <iomanip>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

class ProgressBar {

    private:
        std::ostream & stream;
        uint64_t const total;
        uint64_t processed;
        uint64_t stepWidth;
        unsigned printedSteps;
        std::string prefix;
#ifdef __WITHCDR3THREADS__
        std::mutex MUTEX_updateAndPrintProgress;
#endif

    public:
        void updateAndPrint(uint64_t );
        void print_progress();
        void clear();
        ProgressBar(std::ostream &, uint64_t, unsigned, std::string);
};


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


#endif  // #ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_PROGRESS_BAR_H_
