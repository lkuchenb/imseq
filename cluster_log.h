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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_CLUSTER_LOG_H_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_CLUSTER_LOG_H_

#include <set>
#include <list>
#include <ostream>
#include <iostream>
#include <seqan/file.h>
#include "clone.h"
#include "cdr_utils.h"
#include "segment_meta.h"
#include "cdr3_cli.h"
#include "cluster_result.h"
#include "clone_store.h"

using namespace seqan;

//! A class representing a cluster log entry.
struct ClusterLogEntry {
    //! The larger clone, i.e. the target clone
    Clone<Dna5>         majorClone;
    //! The smaller clone, i.e. the source clone
    Clone<Dna5>         minorClone;
    //! The result record associated with the major clone
    ClusterResult       majorResult;
    //! The result record associated with the minor clone
    ClusterResult       minorResult;
    //! The fraction of the minor clones that are reassigned to the major clones
    double              assignFactor;

    //! Constructor requiring all members to be set
    /*!
    @param      _majorClone     The larger clone, i.e. the target clone
    @param      _majorClone     The larger clone, i.e. the target clone
    @param      _minorResult    The result record associated with the major clone
    @param      _minorResult    The result record associated with the minor clone
    @param      _assignFactor   The fraction of the minor clones that are reassigned to the major clones
    */
    ClusterLogEntry(Clone<Dna5> _majorClone, Clone<Dna5> _minorClone, ClusterResult _majorResult, ClusterResult _minorResult, double _assignFactor)
        : majorClone(_majorClone), minorClone(_minorClone), majorResult(_majorResult), minorResult(_minorResult), assignFactor(_assignFactor) {}
};

inline bool operator<(std::list<ClusterLogEntry>::iterator const & left, std::list<ClusterLogEntry>::iterator const & right) {
    return &(*left) < &(*right);
}

//! A class representing the log of cluster events.
class ClusterLog {

    // ============================================================================
    // TYPEDEFS
    // ============================================================================

    typedef std::list<ClusterLogEntry>::iterator        LogIt;
    typedef std::map<Clone<Dna5>, std::set<LogIt> >     CloneLogEntryMap;

    // ============================================================================
    // DATA MEMBERS
    // ============================================================================

    //! A list of log entries
    std::list<ClusterLogEntry>                          logEntries;

    //! A map for accessing log entries by minor clones
    CloneLogEntryMap                                    logEntriesByMinor;

    //! A map for accessing log entries by major clones
    CloneLogEntryMap                                    logEntriesByMajor;

    // ============================================================================
    // METHODS
    // ============================================================================
    public:

    void logClusterEvent(Clone<Dna5> const & minorClone, String<Clone<Dna5> > const & majorClones, Dna5CloneStore const & cloneStore, String<double> const & reassignedFractions);

    LogIt findEvent(Clone<Dna5> const & minorClone, Clone<Dna5> const & majorClone);

    void validateRates();

    template<typename TGlobal>
    void printClusterLog(TGlobal & global);
    
};

#endif
