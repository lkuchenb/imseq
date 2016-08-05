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

#ifndef IGS_SEGMENT_AMBIGUITY_H
#define IGS_SEGMENT_AMBIGUITY_H

#include <algorithm>
#include "clone.h"
#include "clone_store.h"
#include "cluster_result.h"


template<typename TAlph>
struct ClonesByCDR {
    typedef std::map<String<TAlph>, String<Clone<TAlph> > >     Type;
};

template<typename TClonesByCDR, typename TCloneStore>
void buildClonesByCDR(TClonesByCDR & clonesByCDR, TCloneStore & cloneStore)
{
    clonesByCDR.clear();
    for (typename TCloneStore::iterator it = cloneStore.begin(); it!=cloneStore.end(); ++it)
    {
	appendValue(clonesByCDR[it->first.cdrSeq], it->first);
    }
}

template<typename TSTLContainer>
bool isSubsetOf(TSTLContainer const & major, TSTLContainer const & minor)
{
    if (minor.size() >= major.size())
        return false;
    for (typename TSTLContainer::const_iterator it = minor.begin(); it!=minor.end(); ++it)
        if (major.find(*it) == major.end())
            return false;
    return true;
}

template<typename TSTLContainer>
bool isSubsetOfOrEqual(TSTLContainer const & major, TSTLContainer const & minor)
{
    if (minor.size() > major.size())
        return false;
    for (typename TSTLContainer::const_iterator it = minor.begin(); it!=minor.end(); ++it)
        if (major.find(*it) == major.end())
            return false;
    return true;
}

struct VJReassignment
{
    String<unsigned> minorIDs;
    unsigned majorID, majorVSize, majorJSize;
    VJReassignment() : majorID(0), majorVSize(0), majorJSize(0) {};
    VJReassignment(String<unsigned> _minorIDs, unsigned _majorID, unsigned _majorVSize, unsigned _majorJSize)
	: minorIDs(_minorIDs), majorID(_majorID), majorVSize(_majorVSize), majorJSize(_majorJSize) {};
};

struct 
{
    bool operator()(VJReassignment const & x, VJReassignment const & y)
    {
        return (x.majorVSize + x.majorJSize) > (y.majorVSize + y.majorJSize);
    }
} SortByMinSegSizeFunctor;

/**
 * Compute the relative frequencies from a vector of counts
 */
String<double> computeRelFreqs(String<unsigned> const & counts)
{
    unsigned sum = 0;
    for (Iterator<String<unsigned> const, Rooted>::Type it = begin(counts); !atEnd(it); goNext(it))
	sum += *it;
    String<double> res;
    reserve(res, length(counts));
    for (Iterator<String<unsigned> const, Rooted>::Type it = begin(counts); !atEnd(it); goNext(it))
	 appendValue(res, 1.0 * (*it) / sum);
    return res;
}

template<typename TContainer>
CharString toSeparatedStringX(TContainer const & original, char beginChar='{', char endChar = '}', char sepChar=',') {
    CharString res;
    unsigned reqStrLen = length(original) == 0 ? 0 : length(original) * 2 - 1;
    reqStrLen += 2;
    reserve(res, reqStrLen);

    appendValue(res, beginChar);
    if (length(original) > 0)
        for (typename Iterator<TContainer const>::Type it = begin(original); ; ) {
            appendValue(res, *it);
            ++it;
            if (it==end(original))
                break;
            else
                appendValue(res, sepChar);
        }
    appendValue(res, endChar);

    return(res);
}

/**
 * Given a source count and a vector of relative frequencies, computes the number of
 * clones to be distributed to each target clonotype
 */
String<unsigned> computeTargetCounts(unsigned sourceCount, String<double> const & relFreqs)
{
    String<unsigned> res;
    resize(res, length(relFreqs), 0);
    unsigned redist = 0;
    for (Iterator<String<double> const,Rooted>::Type it = begin(relFreqs); !atEnd(it); goNext(it))
    {
	res[position(it)] = static_cast<unsigned>(std::floor(sourceCount * (*it)));
	redist += res[position(it)];
    }
    redist = sourceCount - redist;
    SEQAN_CHECK(redist <= length(relFreqs), "ERROR [SEQAMB:1001] Please report this error\n");
    for (unsigned x = 0; x < redist; ++x)
	res[x]++;
    return res;
}

template<typename TAlph>
void computeVJReassignmentMap(
	std::map<unsigned, VJReassignment> & reassignmentMap,	// [OUT] the resassignment map
	String<Clone<TAlph> > const & clones)			// [IN ] The clones
{
    reassignmentMap.clear();
    for (unsigned majorID = 0; majorID < length(clones); ++majorID)
    {
        for (unsigned minorID = 0; minorID < length(clones); ++minorID)
        {
            if (majorID == minorID)
                continue;
            if (isSubsetOfOrEqual(clones[majorID].VIds, clones[minorID].VIds) && isSubsetOfOrEqual(clones[majorID].JIds, clones[minorID].JIds)) 
            {
		std::map<unsigned, VJReassignment>::iterator reAssIt = reassignmentMap.find(majorID);
		if (reAssIt != reassignmentMap.end())
		    appendValue(reAssIt->second.minorIDs, minorID);
		else {
		    String<unsigned> minors;
		    appendValue(minors, minorID);
		    VJReassignment r(minors, majorID, (unsigned)clones[majorID].VIds.size(), (unsigned)clones[majorID].JIds.size());
		    reassignmentMap[majorID] = r;
		}
            }
        }
    }
}

template<typename TClone, typename TCloneStore>
void resolveSingleSegAmbiguity(String<TClone> const & clones, TCloneStore & cloneStore)
{
    typedef std::map<unsigned, VJReassignment>	TMap;

    // ===============================================================================
    // Compute the reassignment map, denoting for each clonotype the target clonotypes
    // that it will be redistributed to
    // ===============================================================================

    TMap reassignmentMap;
    computeVJReassignmentMap(reassignmentMap, clones);

    // ===============================================================================
    // Copy the reassignments into a seqan::String and sort them by the total segment
    // count of the major clonotype in descending order. Iterating the reassignments
    // that way ensure that chained reassignments happen in the right order.
    // ===============================================================================

    String<VJReassignment> reassignments;
    for (TMap::const_iterator it = reassignmentMap.begin(); it!=reassignmentMap.end(); ++it)
	appendValue(reassignments, it->second);
    std::sort(begin(reassignments), end(reassignments), SortByMinSegSizeFunctor);

    // ===============================================================================
    // Iterate through the reassignments and execute them
    // ===============================================================================

    for (Iterator<String<VJReassignment> const,Rooted>::Type reAssIt = begin(reassignments); !atEnd(reAssIt); goNext(reAssIt))
    {
	String<unsigned> const & targetIDs = reAssIt->minorIDs;
	// Find all the target clones in the clonestore
	String<typename TCloneStore::iterator> targetCSEIts;
	for (Iterator<String<unsigned> const,Rooted>::Type it = begin(targetIDs); !atEnd(it); goNext(it))
	{
	    appendValue(targetCSEIts, cloneStore.find(clones[*it]));
	    SEQAN_CHECK(targetCSEIts[length(targetCSEIts)-1] != cloneStore.end(), "ERROR [SEQAMB:1004] Please report this error");
	}
	String<unsigned> targetCounts;
	reserve(targetCounts, length(targetIDs));
	for (unsigned minorIndex = 0; minorIndex < length(targetIDs); ++minorIndex)
	    appendValue(targetCounts, targetCSEIts[minorIndex]->second.count);
	// Find the source clone in the clonestore
	typename TCloneStore::iterator sourceCTit = cloneStore.find(clones[reAssIt->majorID]);
	SEQAN_CHECK(sourceCTit!=cloneStore.end(), "ERROR [SEQAMB:1005] Please report this error");
	unsigned sourceCount = sourceCTit->second.count;
	// Compute the redistribution counts
	String<double> targetRelFreqs = computeRelFreqs(targetCounts);
	String<unsigned> redistCounts = computeTargetCounts(sourceCount, targetRelFreqs);
	// Perform the redistribution
	for (unsigned minorIndex = 0; minorIndex < length(targetIDs); ++minorIndex)
	{
	    ClusterResult mergeCR = sourceCTit->second;
	    mergeCR.count = redistCounts[minorIndex];
	    mergeWithClusterResult(targetCSEIts[minorIndex]->second, mergeCR);
	}
	// Remove the old clonotype and it's cluster result from the clone store
	cloneStore.erase(sourceCTit);

    }

    ignoreUnusedVariableWarning(cloneStore);
}

template<typename TCloneStore>
void removeAllelInformation(TCloneStore & cloneStore, CdrReferences const & references)
{
    typedef typename TCloneStore::key_type      TClone;
    TCloneStore newStore;
    for (typename TCloneStore::const_iterator it = cloneStore.begin(); it!=cloneStore.end(); ++it)
    {
	TClone newClone;
	TClone const & oldClone = it->first;
	// Reduce the Clone to primary allel
	for (std::set<unsigned>::const_iterator vSetIt = oldClone.VIds.begin(); vSetIt != oldClone.VIds.end(); ++vSetIt)
	    newClone.VIds.insert(references.leftToFirstAllel[*vSetIt]);
	for (std::set<unsigned>::const_iterator jSetIt = oldClone.JIds.begin(); jSetIt != oldClone.JIds.end(); ++jSetIt)
	    newClone.JIds.insert(references.rightToFirstAllel[*jSetIt]);
	newClone.cdrSeq = oldClone.cdrSeq;
	// Check if the reduced clonotype already exists in the clone store
	typename TCloneStore::iterator newCloneIt = newStore.find(newClone);
	if (newCloneIt != newStore.end())
	    mergeWithClusterResult(newCloneIt->second, it->second);
	else
	    newStore[newClone] = it->second;
    }
    cloneStore = newStore;
}

/**
 * Runs segment ambiguity resolution on the passed clone store
 */
template<typename TCloneStore>
void resolveSegmentAmbiguity(TCloneStore & cloneStore) 
{
    typedef typename TCloneStore::key_type      TClone;
    typedef typename TClone::TAlphabet          TAlph;
    typedef typename ClonesByCDR<TAlph>::Type   TClonesByCDR;

    // Build collection with clones accessible by CDR sequence
    TClonesByCDR clonesByCDR;
    buildClonesByCDR(clonesByCDR, cloneStore);

    std::cerr << "  |-- " << clonesByCDR.size() << " unique CDR3 sequences found.\n";
    

    // Iterate over all clonotypes
    unsigned oldClones = cloneStore.size();
    for (typename TClonesByCDR::iterator it = clonesByCDR.begin(); it!=clonesByCDR.end(); ++it)
    {
        if (length(it->second) < 2)
            continue;
        resolveSingleSegAmbiguity(it->second, cloneStore);
    }
    std::cerr << "  |-- V/J ambiguity resolution reduced from " << oldClones << " to " << cloneStore.size() << " clonotypes.\n";
}

#endif // #ifndef IGS_SEGMENT_AMBIGUITY_H

