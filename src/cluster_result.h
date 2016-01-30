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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_CLUSTER_RESULT_H_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_CLUSTER_RESULT_H_

#include <seqan/sequence.h>
#include "fixed_size_types.h"

using namespace seqan;


class ClusterResult {

    public:
        unsigned count;
        String<uint64_t> sumOfQScores;
        double qMean;
        double qSD;
        ClusterResult();
        String<double> const & getAverageQScores() const;
        String<double> & getWritableAverageQScores();
        void buildAverageQScores();
    private:
        String<double> averageQScores;
        bool wasChanged;
        unsigned lastCount;
};


/*
 * Modifies a base cluster by adding the count of a second cluster as well as adjusting the
 * average qualities accordingly. The average qualities remain on base cluster value at the skip
 * positions
 */
inline void mergeWithClusterResult(ClusterResult & base, ClusterResult const & add, std::set<unsigned> const & skipPos) {
    if (length(base.getWritableAverageQScores()) == 0)
        resize(base.getWritableAverageQScores(), length(add.getAverageQScores()));
    SEQAN_CHECK(length(base.getWritableAverageQScores()) == length(add.getAverageQScores()), 
	    "[ERR] mergeWithClusterResult() encountered ClusterResult instances with different qscore vector lengths. Please report this error.");
    // =1= Adjust the average quality values. Leave out the erroneous positions
    //     in the added result
    std::set<unsigned>::const_iterator skipIt = skipPos.begin();
    for (unsigned i=0; i<length(base.getWritableAverageQScores()); ++i) {
        if (skipIt != skipPos.end() && i==*skipIt) {
            ++skipIt;
            continue;
        }
        base.getWritableAverageQScores()[i] = base.count * base.getWritableAverageQScores()[i] + add.count * add.getAverageQScores()[i];
        base.getWritableAverageQScores()[i] /= base.count + add.count;
    }
    // =2= Increase the count of the base result
    base.count = base.count + add.count;
    // =3= Invalidate the qMean and qSD (will have to be recomputed)
    base.qMean = -1;
    base.qSD = -1;
    // The barcode assignments of the added result are discarded, the ones of
    // the base result remain unchanged
}

inline void mergeWithClusterResult(ClusterResult & base, ClusterResult const & add) {
    mergeWithClusterResult(base, add, std::set<unsigned>());
}

inline void mergeWithClusterResult(ClusterResult & base, ClusterResult const & add, String<unsigned> const & posStr) {
    std::set<unsigned> posSet;
    for (Iterator<String<unsigned> const, Rooted>::Type it = begin(posStr); !atEnd(it); goNext(it))
        posSet.insert(*it);
    mergeWithClusterResult(base, add, posSet);
}

#endif
