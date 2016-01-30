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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_LOGGING_H_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_LOGGING_H_

#include "thread_check.h"
#include <fstream>
#include <iostream>
#include <string>

#ifdef __WITHCDR3THREADS__
#include <mutex>
#endif

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

class Log {
    public:
        typedef std::ostream TStream;
    private:
        std::ofstream * stream;
        std::string path;
        bool hadError;
#ifdef __WITHCDR3THREADS__
        std::mutex lock;
#endif
    public:
        Log(std::string);
        Log();
        ~Log();
        TStream & getStream();
        bool good();
#ifdef __WITHCDR3THREADS__
        void releaseStreamLock();
        void getStreamLock();
#endif
        std::string const & getPath() const;
        void setPath(std::string);
        inline bool wasSpecified() const;
};

inline bool Log::wasSpecified() const {
    return path.length() != 0;
}


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


#endif  // #ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_PROGRESS_BAR_H_
