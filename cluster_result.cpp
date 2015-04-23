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

#include "cluster_result.h"

ClusterResult::ClusterResult() :  count(0), qMean(-1), qSD(-1), wasChanged(false), lastCount(0) {}

void ClusterResult::buildAverageQScores() {
    if (wasChanged) {
        std::cerr << "\n[ERR] Attempt to rebuild averageQScores after they were modified\n";
        exit(1);
    }
    if (count != lastCount) {
        resize(averageQScores, length(sumOfQScores));
        for (unsigned i=0; i<length(sumOfQScores); ++i) 
            averageQScores[i] = 1.0 * sumOfQScores[i] / count;
    }
}

String<double> const & ClusterResult::getAverageQScores() const {
    return averageQScores;
}

String<double> & ClusterResult::getWritableAverageQScores() {
    wasChanged = true;
    return averageQScores;
}

