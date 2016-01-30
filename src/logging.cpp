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

#include "logging.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

Log::Log(std::string p) : stream(NULL), path(p), hadError(false) {}

Log::Log() : stream(NULL), path(std::string("")), hadError(false) {}

Log::~Log() {
    if (stream!=NULL) {
	stream->close();
	delete stream;
    }
}

Log::TStream & Log::getStream() {
    if (stream==NULL)
	stream = new std::ofstream(path.c_str());
    return *stream;
}

bool Log::good() {
    if (path=="" || hadError)
	return false;
    if (stream==NULL)
	stream = new std::ofstream(path.c_str());
    if (!stream->good()) {
	if (!hadError) {
	    hadError = true;
	    std::cerr << "  |xx ERROR - CANNOT WRITE TO FILE '" << path << "'\n";
	}
	return false;
    }
    return true;
}

#ifdef __WITHCDR3THREADS__
void Log::getStreamLock() {
    lock.lock();
}

void Log::releaseStreamLock() {
    lock.unlock();
}
#endif

std::string const & Log::getPath() const {
    return path;
}

void Log::setPath(std::string s) {
    path = s;
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

