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

#include "progress_bar.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

ProgressBar::ProgressBar(std::ostream & s, uint64_t t, unsigned resolution = 100, std::string p = "") :
    stream(s), total(t), processed(0), printedSteps(0), prefix(p) 
{
    stepWidth = total / resolution; 	// Every how many steps do we print
    if (stepWidth == 0) stepWidth = 1;
}	


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


/**
 * Adds the finished count to the processed count and re-prints the
 * progressbar if the threshold for re-printing is reached.
 * THREADSAFE
 */
void ProgressBar::updateAndPrint(uint64_t finished) {
#ifdef __WITHCDR3THREADS__
    std::unique_lock<std::mutex> lock(this->MUTEX_updateAndPrintProgress);
#endif
    processed += finished;
    if (processed > total)
        processed = total;
    unsigned newSteps = processed / stepWidth;
    if (newSteps > printedSteps) {
	this->printedSteps = newSteps;
	this->print_progress();
    }

}

void ProgressBar::print_progress() {
    double p = 100.0 * processed / total;
    stream << prefix << "[";
    for (unsigned i=1; i<=50; ++i) {
        if (i==24) {
            stream << std::setprecision(2) << std::fixed << p;
            i+= p<10 ? 3 : 4;
            continue;
        }
        if (i<p/2)
            stream << "#";
        else
            stream << " ";
    }
    stream << "]\r";
}

void ProgressBar::clear() {
    std::string prefixWhitespace;
    prefixWhitespace.resize(prefix.size(), ' ');
    stream << prefixWhitespace;
    stream << "                                                        \r";
}

