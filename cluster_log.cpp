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

#include "cluster_log.h"

ClusterLog::LogIt ClusterLog::findEvent(Clone<Dna5> const & minorClone, Clone<Dna5> const & majorClone) {
    std::set<LogIt> minorEvents = logEntriesByMinor[minorClone];
    std::set<LogIt> majorEvents = logEntriesByMajor[majorClone];
    std::set<LogIt> overlap = setIntersect(minorEvents, majorEvents);
    SEQAN_CHECK(overlap.size() == 1 || overlap.empty(), "ERROR 1005 - findEvent() redundancy. Please report this error.");
    if (overlap.empty())
	return logEntries.end();
    else
	return *(overlap.begin());
}

void ClusterLog::logClusterEvent(Clone<Dna5> const & minorClone, String<Clone<Dna5> > const & majorClones, Dna5CloneStore const & cloneStore, String<double> const & reassignedFractions) {

    SEQAN_CHECK(seqanStringSum(reassignedFractions) - 1 < 0.000000001, "ERROR 1003 - logClusterEvent() consistency check. Please report this error!");
    SEQAN_CHECK(length(majorClones) == length(reassignedFractions), "ERROR 1004 - logClusterEvent() size inconsistency. Please report this error!");
    SEQAN_CHECK(logEntriesByMinor.find(minorClone)==logEntriesByMinor.end(), "ERROR 1005 - logClusterEvent() clustering a minor twice. Please report this error!");

    // Look up the major clones results
    std::vector<ClusterResult const *> majorResultPntrs;
    reserve(majorResultPntrs, length(majorClones));
    for (Iterator<String<Clone<Dna5> > const,Rooted>::Type it = begin(majorClones); !atEnd(it); goNext(it))
	majorResultPntrs.push_back(&(mapMustHaveKey(cloneStore, *it, "ERROR 1009 - logClusterEvent()", 1009)));

    // Check if the minor clone was already a major clone before
    // and reassign the pre-minors
    std::set<LogIt> toMinorEvents = logEntriesByMajor[minorClone];
    for (std::set<LogIt>::iterator logEntryItIt = toMinorEvents.begin(); logEntryItIt != toMinorEvents.end(); ++logEntryItIt) {
	LogIt logEntryIt = *logEntryItIt;
	for (unsigned majId = 0; majId<length(majorClones); ++majId) {
	    LogIt oldMinorNewMajor = findEvent(logEntryIt->minorClone, majorClones[majId]);
	    // Case 1: Pre-minor already has a fraction assigned to this major
	    if (oldMinorNewMajor != logEntries.end()) {
		oldMinorNewMajor->assignFactor += reassignedFractions[majId] * logEntryIt->assignFactor;
	    }
	    // Case 2: A new pre-minor <-> major event has to be created
	    else {
		// Find the result associated with the major clone
		logEntries.push_front(ClusterLogEntry(majorClones[majId], logEntryIt->minorClone, *majorResultPntrs[majId], logEntryIt->minorResult, reassignedFractions[majId] * logEntryIt->assignFactor));
		ClusterLog::LogIt clePnt = logEntries.begin();
		logEntriesByMinor[clePnt->minorClone].insert(clePnt);
		logEntriesByMajor[clePnt->majorClone].insert(clePnt);
	    }
	}
	// Erase the intermediate event
	size_t x = logEntriesByMinor[logEntryIt->minorClone].erase(logEntryIt);	// Erase by key_type
	size_t y = logEntriesByMajor[logEntryIt->majorClone].erase(logEntryIt);	// Erase by key_type
	if (x != 1 || y != 1) {
	    std::cerr << "\n\nERROR 1002 - logClusterEvent() [" << x << ";" << y << "] - Please report this error!\n\n";
	    exit(1002);
	}
	logEntries.erase(logEntryIt);						// Erase by iterator
    }
    ClusterResult const & minorResult = mapMustHaveKey(cloneStore, minorClone, "ERROR 1010 - logClusterEvent()", 1010);
    // Create the actual log events
    for (unsigned majId = 0; majId<length(majorClones); ++majId) {
	logEntries.push_front(ClusterLogEntry(majorClones[majId], minorClone, *majorResultPntrs[majId], minorResult, reassignedFractions[majId]));
	ClusterLog::LogIt clePnt = logEntries.begin();
	logEntriesByMinor[minorClone].insert(clePnt);
	logEntriesByMajor[majorClones[majId]].insert(clePnt);
    }
}

void ClusterLog::validateRates() {
    for (CloneLogEntryMap::const_iterator minorEntriesIt = logEntriesByMinor.begin(); minorEntriesIt!=logEntriesByMinor.end(); ++minorEntriesIt) {
	double sum = 0;
	std::stringstream ss;
	for (std::set<LogIt>::const_iterator logItIt = minorEntriesIt->second.begin(); logItIt != minorEntriesIt->second.end(); ++logItIt) {
	    sum += (*logItIt)->assignFactor;
	    ss << (*logItIt)->assignFactor << "  ";
	}
	std::cout << "ASSIGN FACTOR SUM " << sum << "=" << ss.str() << std::endl;
    }
}



