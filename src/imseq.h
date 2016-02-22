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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_CDR3FINDER_H_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_CDR3FINDER_H_

#include <cstring>
#include <cstdlib>
#include <string>
#include <sstream>
#include <utility>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_align.h>

#include "aa_translate.h"
#include "progress_bar.h"
#include "thread_check.h"
#include "cluster_result.h"
#include "fixed_size_types.h"
#include "collection_utils.h"
#include "segment_meta.h"
#include "cdr_utils.h"
#include "clone.h"
#include "cluster_log.h"
#include "vjMatching.h"
#include "overlap_specs.h"
#include "referencePreparation.h"
#include "timeFormat.h"
#include "file_utils.h"
#include "sequence_data.h"
#include "segment_ambiguity.h"
#include "fastq_io.h"
#include "fastq_multi_record.h"
#include "reject.h"
#include "barcode_correction.h"
#include "input_information.h"

#ifdef __WITHCDR3THREADS__
#include <mutex>
#include <thread>
#include "thread_pool.h"
#include "reject.h"
#endif

#ifndef REJECTLOG
#define REJECTLOG if (false) std::cerr
#endif

using namespace seqan;

#ifdef FULLREADOUT
StringSet<CharString> __READIDS;
StringSet<CharString> __VGENES;
StringSet<CharString> __JGENES;
StringSet<String<AminoAcid> > __CDRSEQS;
StringSet<String<Dna5Q> > __FULLSEQS;
#endif

void PressEnterToContinue()
{
    std::cerr << "Press ENTER to continue... " << std::flush;
    std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );
}


// ============================================================================
// Locks
// ============================================================================

#ifdef __WITHCDR3THREADS__
std::mutex MUTEX_ofstream_rejectlog;
std::mutex STDCOUT_LOCK;
#endif

// ============================================================================
// Forwards
// ============================================================================

template<typename TRow>
inline unsigned viewEndPosition(TRow const &);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


struct AnalysisResult{
    Clone<Dna5>     clone;
    RejectReason    reject;
    CharString      fullOutSuffix;
    String<uint64_t>     cdrQualities;
    
    AnalysisResult(Clone<Dna5> _clone, CharString _fullOutSuffix, String<uint64_t> _cdrQualities) : clone(_clone), fullOutSuffix(_fullOutSuffix), cdrQualities(_cdrQualities) { 
        reject = NONE;
    }

    AnalysisResult(RejectReason r) : reject(r) { }

    AnalysisResult() {
        reject = NONE;
    }
};

inline bool operator==(const AnalysisResult& lhs, const AnalysisResult& rhs)
{
    //std::cerr << (lhs.clone==rhs.clone) << " " << (lhs.reject == rhs.reject) << " " <<  (lhs.fullOutSuffix == rhs.fullOutSuffix) << std::endl;
    bool res = (lhs.clone==rhs.clone && lhs.reject == rhs.reject && lhs.fullOutSuffix == rhs.fullOutSuffix);
    if (!res) {
        std::cerr << "\n-----------\nFULLOUTSUFFIX\n";
        std::cerr << lhs.fullOutSuffix << " == " << rhs.fullOutSuffix << " = " << (lhs.fullOutSuffix== rhs.fullOutSuffix) << std::endl;
        std::cerr << "\n-----------\nREJECT\n";
        std::cerr << _CDRREJECTS[lhs.reject] << " == " << _CDRREJECTS[rhs.reject] << " = " << (lhs.reject== rhs.reject) << std::endl;
        std::cerr << "\n-----------\nCLONE\n";
        std::cerr << toString(lhs.clone) << " == " << toString(rhs.clone) << " = " << (lhs.clone==rhs.clone) << std::endl;
    }
    return res;
}

inline bool operator!=(const AnalysisResult& lhs, const AnalysisResult& rhs)
{
    return (!(lhs==rhs));
}

struct CacheEntry {
    AnalysisResult analysisResult;
    unsigned cdrBeginPos;
    unsigned cdrEndPos;
    String<unsigned> modifiedPositions;

    CacheEntry(AnalysisResult const & ar, unsigned _cdrBeginPos = 0, unsigned _cdrEndPos = 0, String<unsigned> const _modifiedPositions = String<unsigned>()) : cdrBeginPos(_cdrBeginPos), cdrEndPos(_cdrEndPos), modifiedPositions(_modifiedPositions) {
        analysisResult = ar;
        if (ar.reject!=NONE) {
            cdrBeginPos = 0;
            cdrEndPos = 0;
            clear(modifiedPositions);
        }
        // As the cache entries are *sequence* specific, but not *read* specific, the quality of the
        // cdr3 region is not cached. It has to be recomputed for each read individually.
        resize(analysisResult.cdrQualities, 0);
    }

};

inline bool operator==(const CacheEntry& lhs, const CacheEntry& rhs)
{
    bool res = (lhs.analysisResult == rhs.analysisResult && lhs.cdrBeginPos == rhs.cdrBeginPos && lhs.cdrEndPos == rhs.cdrEndPos);
    if (!res) {
        std::cerr << "\nlhs.cdrBeginPos=" << lhs.cdrBeginPos
            << "; rhs.cdrBeginPos=" << rhs.cdrBeginPos
            << "; lhs.cdrEndPos=" << lhs.cdrEndPos
            << "; rhs.cdrEndPos=" << rhs.cdrEndPos << "\n";
    }
    return(res);
}

inline bool operator!=(const CacheEntry& lhs, const CacheEntry& rhs)
{
    return (!(lhs==rhs));
}


typedef std::map<String<Dna5>, CacheEntry>    TCacheStore;
#ifdef __WITHCDR3THREADS__
std::mutex CACHE_STORE_MUTEX;
#endif

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Typedefs
// ============================================================================

typedef std::map<Clone<Dna5>, ClusterResult>    TCloneStore;

/**
 * Transforms a String<T> into a String<char>, making use of convert() and
 * using specified characters for element separation as well as begin / end
 * notation
 */
template<typename TContainer>
CharString toSeparatedString(TContainer const & original, char beginChar='{', char endChar = '}', char sepChar=',') {
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


template<typename TContainer>
std::string toSeparatedStringConverted(TContainer const & original, char beginChar='{', char endChar = '}', char sepChar=',') {
    std::stringstream ss;

    ss << beginChar;
    if (length(original) > 0)
        for (typename Iterator<TContainer const>::Type it = begin(original); ; ) {
            ss << *it;
            ++it;
            if (it==end(original))
                break;
            else
                ss << sepChar;
        }
    ss << endChar;

    return(ss.str());
}


inline std::string toStringA(AnalysisResult const & ar) {
    std::stringstream ss;
    ss << "AnalysisResult { clone=" << toString(ar.clone) << "; reject=" << _CDRREJECTS[ar.reject] << "; fullOutSuffix=" << "unknown" << "; cdrQualities=" << ar.cdrQualities;
    return(ss.str());
}

inline std::string toString(CacheEntry const & ce) {
    std::stringstream ss;
    ss <<  "CacheEntry { analysisResult=" << toStringA(ce.analysisResult) << "; cdrBeginPos=" << ce.cdrBeginPos << "; cdrEndPos="
        << ce.cdrEndPos << "; modifiedPositions=" << toSeparatedString(ce.modifiedPositions) << "}";
    return ss.str();
}

template<typename T>
std::ostream& operator<<(std::ostream& os, Clone<T> const & clone) {
    os << "[left: " << toSeparatedStringSTL(clone.VIds) << "; seq:" << clone.cdrSeq << "; right: " << toSeparatedStringSTL(clone.JIds) << "]";
    return(os);
}

std::ostream& operator<<(std::ostream& os, ClusterResult const & res) {
    os << "[count:" << res.count << "; avgquals: " << toSeparatedStringConverted(res.getAverageQScores(), '<', '>') <<  "; mean:" << res.qMean << "; sd:" << res.qSD << "]";
    return(os);
}

struct ClonesBySize {

    TCloneStore const & cloneStore;

    bool operator() (const Clone<Dna5> * lhs, const Clone<Dna5> * rhs) const {
        TCloneStore::const_iterator _lhs = cloneStore.find(*lhs);
        TCloneStore::const_iterator _rhs = cloneStore.find(*rhs);
        return _lhs->second.count < _rhs->second.count;
    }

    ClonesBySize(TCloneStore const & cs) : cloneStore(cs) {};
};


// ============================================================================
// Functions
// ============================================================================

/**
 * Ugly wrapper ignoring the specified mode string
 */
void open(std::ifstream & is, const char * filename, CharString const & mode) {
    ignoreUnusedVariableWarning(mode);
    is.open(filename, std::ios_base::in | std::ios_base::binary);
}


template<typename TContainer>
inline std::pair<double,double> meanSd(TContainer const & c) {
    double mean = 0;
    for (typename Iterator<TContainer const, Rooted>::Type it = begin(c); !atEnd(it); goNext(it))
        mean += *it;
    mean /= length(c);

    double sd = 0;
    for (typename Iterator<TContainer const, Rooted>::Type it = begin(c); !atEnd(it); goNext(it)) {
        double x = *it - mean;
        sd += x*x;
    }
    sd = sd / (length(c) - 1);
    sd = std::sqrt(sd);
    return(std::make_pair(mean,sd));
}


void close(std::ifstream & is) {
    is.close();
}

seqan::SeqFileIn & getInfoStream(SeqInputStreams<SingleEnd> & input) {
    return input.stream;
}

seqan::SeqFileIn & getInfoStream(SeqInputStreams<PairedEnd> & input) {
    return input.fwStream;
}

void resetStreams(SeqInputStreams<SingleEnd> & input) {
    close(input.stream);
    open(input.stream, input.path.c_str());
}

void resetStreams(SeqInputStreams<PairedEnd> & input) {
    close(input.fwStream);
    close(input.revStream);
    open(input.fwStream, input.fwPath.c_str());
    open(input.revStream, input.revPath.c_str());
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// CLASS: ClusterPair
////////////////////////////////////////////////////////////////////////////////////////////////////

class ClusterPair {

    public:

        Clone<Dna5> const * minorPnt;
        Clone<Dna5> const * majorPnt;
        String<unsigned> minorErrPositions;
        bool isValidVal;
        std::string reason;

    public:

        ClusterPair(std::pair<Clone<Dna5> const,ClusterResult> const &, std::pair<Clone<Dna5> const, ClusterResult> const &, unsigned maxErrors, std::string const & );
        bool isValid();
        Clone<Dna5> const & getMajorClone() const;
        Clone<Dna5> const & getMinorClone() const;
        String<unsigned> const & getMinorErrPositions() const;
};


ClusterPair::ClusterPair(TCloneStore::value_type const & p1,
        TCloneStore::value_type const & p2,
        unsigned maxErrors = -1u, std::string const & _reason = "") : isValidVal(true), reason(_reason)
{
    if (p1.second.count > p2.second.count) {
        this->majorPnt = &(p1.first);
        this->minorPnt = &(p2.first);
    } else {
        this->majorPnt = &(p2.first);
        this->minorPnt = &(p1.first);
    }
   
    if (length(minorPnt->cdrSeq) != length(majorPnt->cdrSeq)) {
        std::cerr << "[ERR] Attempt to create a ClusterPair object from clusters with different CDR3 lengths. Please report this error." << std::endl;
        exit(1);
    }

    for (unsigned i=0; i<length(minorPnt->cdrSeq); ++i) {
        if (minorPnt->cdrSeq[i] != majorPnt->cdrSeq[i]) {
            appendValue(minorErrPositions, i);
            if (length(minorErrPositions) > maxErrors) {
                isValidVal = false;
                break;
            }
        } 
    }
}

bool ClusterPair::isValid() {
    return this->isValidVal;
}

Clone<Dna5> const & ClusterPair::getMajorClone() const {
    return *(this->majorPnt);
}

Clone<Dna5> const & ClusterPair::getMinorClone() const {
//    std::cerr << "getMinorClone() called! Returned: " << *(this->minorPnt) << '\n';
    return *(this->minorPnt);
}

String<unsigned> const & ClusterPair::getMinorErrPositions() const {
    return this->minorErrPositions;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
inline bool sharesMember(std::set<T> const & set1, std::set<T> const & set2) {
    for (typename std::set<T>::const_iterator it1 = set1.begin(); it1!=set1.end(); ++it1)
        for (typename std::set<T>::const_iterator it2 = set2.begin(); it2!=set2.end(); ++it2)
            if (*it1 == *it2) return true;
    return false;
}

/**
 * Checks if all the error positions in a cluster pair correspond to sufficiently
 * low quality values in the minor cluster
 */
inline bool checkSDDiff(CdrOptions const & options, ClusterPair const & cp, TCloneStore const & cloneStore) {
    // Find the minor cluster result
    TCloneStore::const_iterator cse = cloneStore.find(cp.getMinorClone());

    if (cse == cloneStore.end())
        std::cerr << "[ERR] checkSDDiff() failed to find a clone. Please report this error." << std::endl;
    else if (cse->second.qMean == -1 || cse->second.qSD == -1)
        std::cerr << "[ERR] checkSDDiff() encountered a cluster result without mean / sd data. Please report this error." << std::endl;

    double maxValue = cse->second.qMean - (1.0 * options.minSdDevi * cse->second.qSD);
    //std::cerr << "[maxValue " << maxValue << "] ";
    for (Iterator<String<unsigned> const,Rooted>::Type errPos = begin(cp.getMinorErrPositions()); !atEnd(errPos); goNext(errPos)) {
        //std::cerr << cse->second.getAverageQScores()[*errPos] << " ";
        if (cse->second.getAverageQScores()[*errPos] > maxValue)
            return false;
    }
    //std::cerr << std::endl;

    return true;
}

template<typename T>
void addAllElements(std::set<T> & target, std::set<T> const & source) {
    for (typename std::set<T>::const_iterator it = source.begin(); it!=source.end(); ++it)
        target.insert(*it);
}


/**
 * Determines disjoint groups from sets, i.e. puts those sets together in one group that
 * are connected by non-empty intersection, possibly with intermediate groups
 */
void buildGroups(String<std::set<unsigned> > & groups, String<std::set<unsigned> > const & matchSet) {
    
    String<std::set<unsigned> > tmpGroups;
    for (Iterator<String<std::set<unsigned> > const,Rooted>::Type matchSetElem = begin(matchSet); !atEnd(matchSetElem); goNext(matchSetElem)) {
        //std::cerr << "BEFORE " << tmpGroups << std::endl;
        //std::cerr << *matchSetElem << " [MatchSetElem]" << std::endl;
        if (length(tmpGroups)==0) {
            resize(tmpGroups,1);
            tmpGroups[0].insert(position(matchSetElem));
            //std::cerr << "AFTER  " << tmpGroups << "\n" << std::endl;
            continue;
        }

        // Find all groups that contain matchsets that intersect with this one
        String<unsigned> shares;
        for (int i = 0; i < static_cast<int>(length(tmpGroups)); ++i) {
            for (std::set<unsigned>::const_iterator gId = tmpGroups[i].begin(); gId!=tmpGroups[i].end(); ++gId) 
                if (sharesMember(matchSet[*gId], *matchSetElem))  {
                    appendValue(shares, i);
                    break;
                }
        }

        if (length(shares)>0) {
            // Merge groups that are connected by the current match-set
            Iterator<String<unsigned>, Rooted>::Type addMe = begin(shares);
            for (goNext(addMe); !atEnd(addMe); goNext(addMe)) {
                addAllElements(tmpGroups[shares[0]], tmpGroups[*addMe]);
                tmpGroups[*addMe].clear();
            }
            tmpGroups[shares[0]].insert(position(matchSetElem));
        } else {
            resize(tmpGroups,length(tmpGroups)+1);
            tmpGroups[length(tmpGroups)-1].insert(position(matchSetElem));
        }
        //std::cerr << "AFTER  " << tmpGroups << "\n" << std::endl;
    }
    clear(groups);
    for (Iterator<String<std::set<unsigned> >,Rooted>::Type it = begin(tmpGroups); !atEnd(it); goNext(it))
        if (it->size() != 0)
            appendValue(groups, *it);
}

template<typename T>
std::set<T> setIntersection(std::set<T> const & setA, std::set<T> const & setB) {
    std::set<T> res;
    for (typename std::set<T>::const_iterator aIt = setA.begin(); aIt!=setA.end(); ++aIt)
        if (setB.find(*aIt) != setB.end())
            res.insert(*aIt);
    return(res);
}

void buildConsensusGroups(String<std::set<unsigned> > & consensus, String<std::set<unsigned> > const & groupA, String<std::set<unsigned> > const & groupB) {
    typedef String<std::set<unsigned> >  TGroup;
    
    clear(consensus);
    for (Iterator<TGroup const, Rooted>::Type aIt = begin(groupA); !atEnd(aIt); goNext(aIt)) {
        for (Iterator<TGroup const, Rooted>::Type bIt = begin(groupB); !atEnd(bIt); goNext(bIt)) {
            std::set<unsigned> inter = setIntersection(*aIt, *bIt);
            if (inter.size()>0)
                appendValue(consensus, inter);
        }
    }

}

template<typename TAlphabet, typename TCloneStore_>
void mergeClones(String<Clone<TAlphabet> > const & targetClones, Clone<TAlphabet> const & newClone, TCloneStore_ & cloneStore) {
    if (length(targetClones)<=1)
        return;
    // Create an "empty" cluster result (count 0)
    ClusterResult cluRes;
    for (typename Iterator<String<Clone<TAlphabet> > const,Rooted>::Type it = begin(targetClones); !atEnd(it); goNext(it)) {
        typename TCloneStore_::iterator oldElem = cloneStore.find(*it);
        if (oldElem == cloneStore.end()) {
            std::cerr << "[ERR] mergeClones() failed to find cluster result. Please report this error!" << std::endl;
            exit(1);
        }
        mergeWithClusterResult(cluRes, oldElem->second);
        //std::cout << "DEL " << oldElem->first << " - " << oldElem->second << std::endl;
        cloneStore.erase(oldElem);
    }
    cloneStore[newClone] = cluRes;
    //std::cout << "ADD " << newClone << " - " << cluRes << "\n\n" << std::flush;
}

template<typename TCloneStore_>
void mergeIdenticalCDRs(TCloneStore_ & cloneStore) { 
    typedef std::map<String<Dna5>, String<typename TCloneStore_::value_type> >   TCloneStoreBySeq;
    typedef std::map<std::set<unsigned>, unsigned>                               TCounterMap;
    typedef typename TCloneStore_::key_type::TAlphabet                           TAlphabet;

    std::cerr << "===== Merging identical CDR3s with ambiguous segment matches" << std::endl;

    // ============================================================================
    // Make all clonestore elements accessible by their sequence
    // ============================================================================

    TCloneStoreBySeq cloneStoreBySeq;
    for (typename TCloneStore_::const_iterator cse = cloneStore.begin(); cse!=cloneStore.end(); ++cse)
        appendValue(cloneStoreBySeq[cse->first.cdrSeq], *cse);

    std::cerr << "  |-- " << cloneStoreBySeq.size() << " unique CDR3 sequences found" << std::endl;

    for (typename TCloneStoreBySeq::const_iterator cbsIt = cloneStoreBySeq.begin(); cbsIt != cloneStoreBySeq.end(); ++cbsIt) {
        if (length(cbsIt->second)<=1)
            continue;

        // ============================================================================
        // Build the sets of left / right match ids
        // ============================================================================

        String<std::set<unsigned> > leftMatchSets, rightMatchSets;
        for (typename Iterator<typename TCloneStoreBySeq::mapped_type const, Rooted>::Type cse = begin(cbsIt->second); !atEnd(cse); goNext(cse)) {
            appendValue(leftMatchSets, cse->first.VIds);
            appendValue(rightMatchSets, cse->first.JIds);
        }

        // ============================================================================
        // Group the clone store elements based on their left / right  matches. The 
        // clone store elements are accessed by their in index in the String<> that
        // contains them within the cloneStoreBySeq data structure.
        // ============================================================================

        String<std::set<unsigned> > leftGroups, rightGroups;
        //std::cerr << "LEFT" << std::endl;
        buildGroups(leftGroups, leftMatchSets);
        //std::cerr << "RIGHT" << std::endl;
        buildGroups(rightGroups, rightMatchSets);

        // ============================================================================
        // Compute the consensus grouping
        // ============================================================================

        String<std::set<unsigned> > consensusGroups;
        buildConsensusGroups(consensusGroups, leftGroups, rightGroups);

        // ============================================================================
        // In each group, use the intersection as the new segment set or the union if
        // empty [edit: two alternative methods are currently evaluated]
        // ============================================================================

        for (Iterator<String<std::set<unsigned> >,Rooted>::Type group = begin(consensusGroups); !atEnd(group); goNext(group)) {
            if ( group->size() <= 1 )
                continue;
            std::set<unsigned> leftConsensus, rightConsensus, leftUnion, rightUnion;
            String<Clone<TAlphabet> > clones;
            TCounterMap leftMatchCounter, rightMatchCounter;

            for (std::set<unsigned>::const_iterator indexIt = group->begin(); indexIt!=group->end(); ++indexIt) {
                appendValue(clones, cbsIt->second[*indexIt].first);
                // Counting... later the top ranked set is chosen
                leftMatchCounter[leftMatchSets[*indexIt]] += cbsIt->second[*indexIt].second.count;
                rightMatchCounter[rightMatchSets[*indexIt]] += cbsIt->second[*indexIt].second.count;
            }

            // Counting method
            unsigned leftMaxCnt = 0, rightMaxCnt = 0;
            // (1) find max
            for (TCounterMap::const_iterator it = leftMatchCounter.begin(); it!=leftMatchCounter.end(); ++it)
                if (it->second > leftMaxCnt) {
                    leftMaxCnt    = it->second;
                    leftConsensus = it->first;
                }
            // (2) Intersect if possible TODO [!!!] This might be erroneous! Result depends on the order!
            for (TCounterMap::const_iterator it = leftMatchCounter.begin(); it!=leftMatchCounter.end(); ++it) {
                std::set<unsigned> intersect = setIntersection(leftConsensus, it->first);
                if (intersect.size() > 0)
                    leftConsensus = intersect;
            }
            // Repeat for right segment...
            for (TCounterMap::const_iterator it = rightMatchCounter.begin(); it!=rightMatchCounter.end(); ++it)
                if (it->second > rightMaxCnt) {
                    rightMaxCnt    = it->second;
                    rightConsensus = it->first;
                }
            for (TCounterMap::const_iterator it = rightMatchCounter.begin(); it!=rightMatchCounter.end(); ++it) {
                std::set<unsigned> intersect = setIntersection(rightConsensus, it->first);
                if (intersect.size() > 0)
                    rightConsensus = intersect;
            }

            // END Counting method

            Clone<TAlphabet> c = { leftConsensus, cbsIt->first, rightConsensus };
            mergeClones(clones, c, cloneStore);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// CLASS: CloneToPairsMap
////////////////////////////////////////////////////////////////////////////////////////////////////

 /******************************************************************************************
 * Wrapper around a std::map required for the quality clustering. The ClusterPair objects
 * for each Clone are stored in decrementing order of the major cluster size.
 ******************************************************************************************/
class CloneToPairsMap {
    private:

        TCloneStore const & cs;

        struct Comp {
            
            TCloneStore const & cs;

            bool operator() (const ClusterPair & lhs, const ClusterPair & rhs) const {
                // This !less implementation must ensure that even if two CRs have the same
                // major count, they will still always be ordered in the same way in a 
                // multiset. Therefore, the comparision is extended to take the entire
                // major clone into account (the minor clone does not need to be taken into
                // account since it is always the same in this store).
                ClusterResult const & _lhs = cs.find(lhs.getMajorClone())->second;
                ClusterResult const & _rhs = cs.find(rhs.getMajorClone())->second;
                return _lhs.count > _rhs.count || ((_lhs.count == _rhs.count) && (lhs.getMajorClone() < rhs.getMajorClone()));
            }

            Comp(TCloneStore const & _cs) : cs(_cs) {};

        };

        std::map<Clone<Dna5>, std::multiset<ClusterPair, Comp> > store;

    public:


        typedef std::map<Clone<Dna5>, std::multiset<ClusterPair, Comp> > TStore;

        CloneToPairsMap(TCloneStore const & cloneStore) : cs(cloneStore) {};

        TStore::mapped_type & operator[] ( const TStore::key_type& x ) {
            if (store.find(x) == store.end())
                store.insert(std::make_pair(x, TStore::mapped_type(Comp(cs))));

            return store.find(x)->second;
        }

        TStore const & getMap() {
            return store;
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////////////////////////

inline unsigned totalNumberOfClones(TCloneStore const & cloneStore) {
    unsigned count = 0;
    for (TCloneStore::const_iterator it = cloneStore.begin(); it != cloneStore.end(); ++it) {
        count += it->second.count;
    }
    return count;
}


#ifdef __WITHCDR3THREADS__
std::mutex MUTEX_getLQPositions;
#endif
typedef std::map<const ClusterResult*, std::set<unsigned> > TLQPositionStore;
inline std::set<unsigned> const getLQPositions(TLQPositionStore & store, ClusterResult const & cluRes, CdrOptions const & options) {

#ifdef __WITHCDR3THREADS__
    std::unique_lock<std::mutex> lock(MUTEX_getLQPositions);
#endif

    TLQPositionStore::const_iterator it = store.find(&cluRes);
    if (it != store.end())
        return(it->second);

    std::pair<double,double> meanSDVals = meanSd(cluRes.getAverageQScores());
    double maxVal = meanSDVals.first - options.minSdDevi * meanSDVals.second;
    std::set<unsigned> positions;
    for (unsigned i=0; i<length(cluRes.getAverageQScores()); ++i)
        if (cluRes.getAverageQScores()[i] <= maxVal)
            positions.insert(i);

    store[&cluRes] = positions;

    return positions;
}

/**
 * Only the major component of the pairs is taken into account
 */
inline void computeRedistCounts(String<unsigned> & redistributionCounts, String<double> & redistributionRatios, Clone<Dna5> const & minorClone, String<ClusterPair> const & targetPairs, TCloneStore const & cloneStore) {
    // Check if there is anything to do
    if (length(targetPairs)==0)
        return; 

    // Clear and resize the target containers
    clear(redistributionCounts);
    clear(redistributionRatios);
    resize(redistributionCounts, length(targetPairs));
    resize(redistributionRatios, length(targetPairs));

    // Compute the total number of reads within the target clonotypes
    unsigned sum = 0;
    for (Iterator<String<ClusterPair> const,Rooted>::Type it = begin(targetPairs); !atEnd(it); goNext(it)) {
        TCloneStore::const_iterator cse = cloneStore.find(it->getMajorClone());
        if (cse == cloneStore.end()) {
            std::cerr << "[ERR] error in computeRedistCounts() [1]" << std::endl;
            exit(1);
        }
        sum += cse->second.count;
    }

    TCloneStore::const_iterator cse = cloneStore.find(minorClone);

    SEQAN_CHECK(cse != cloneStore.end(), "ERROR 1006 - computeRedistCounts(). Please report this error.");

    unsigned minorCount = cse->second.count;
    unsigned usedMinorCount = minorCount;

    // Compute the ratio of the reads within the target
    // clonotypes
    for (Iterator<String<ClusterPair> const,Rooted>::Type it = begin(targetPairs); !atEnd(it); goNext(it)) {
        redistributionRatios[position(it)] = 1.0 * cloneStore.find(it->getMajorClone())->second.count / sum;
        unsigned nClones = redistributionRatios[position(it)] * minorCount;
        redistributionCounts[position(it)] = nClones;
        usedMinorCount -= nClones;
    }

    for (Iterator<String<unsigned>,Rooted>::Type it = begin(redistributionCounts); !atEnd(it); goNext(it)) {
        if (usedMinorCount<=0)
            break;
        ++(*it);
        --usedMinorCount;
    }

    if (usedMinorCount>0) {
        std::cerr << "[ERR] error in computeRedistCounts() [3]" << std::endl;
        exit(1);
    }
}

template<typename TVal, typename TSpec>
double meanVal(String<TVal, TSpec> const & s) {
    double res = 0;
    for (typename Iterator<String<TVal,TSpec> const>::Type it = begin(s); it!=end(s); ++it) 
        res += *it;
    if (length(s) > 0) res /= length(s);
    return(res);
}

void dropLowQualityClusters(TCloneStore & cloneStore, double const minMeanQuality) {


    // Check if there is work to be done
    if (minMeanQuality <= 0)
        return;

    std::cerr << "===== Removal of low quality clonotype clusters" << std::endl;
    std::cerr << "  |-- User defined threshold: " << minMeanQuality << std::endl;

    typedef String<TCloneStore::key_type>   TRemoveString;

    // Store all CTs that violate the quality condition in a container
    TRemoveString toRemove;
    for (TCloneStore::const_iterator cse = cloneStore.begin(); cse!=cloneStore.end(); ++cse) 
        if(meanVal(cse->second.getAverageQScores()) < minMeanQuality)
            appendValue(toRemove, cse->first);

    std::cerr << "  |-- " << length(toRemove) << " out of " << length(cloneStore) << " clonotype clusters affected." << std::endl;

    // Remove them from the clone store
    for (Iterator<TRemoveString>::Type it = begin(toRemove); it!=end(toRemove); ++it) 
        cloneStore.erase(*it);

}

struct ClusterStats {
    uint64_t falseNegatives;
    uint64_t trueNegatives;
    uint64_t checkedPairs;
    uint64_t alignmentCount;
    ClusterStats() : falseNegatives(0), trueNegatives(0), checkedPairs(0), alignmentCount(0) {}
};

enum ClusterFailure {
    NO_FAILURE          = 0,
    VJ_MISMATCH         = 1,
    CDR3LEN_MISMATCH    = 2,
    TOO_MANY_ERRORS     = 3,
    NO_QUAL_DEVIATION   = 4,
    RATIO_EXCEEDED      = 5
};
std::string const CLUSTER_FAIL_REASONS[6] = { "NONE","VJ_MISMATCH","CDR3LEN_MISMATCH","TOO_MANY_ERRORS", "NO_QUAL_DEVIATION", "RATIO_EXCEEDED" };

template<typename T>
unsigned maxLength(T const & a, T const & b) {
    return(std::max(length(a), length(b)));
}

#ifdef __WITHCDR3THREADS__
std::mutex MUTEX_findClusterMates_store_result, MUTEX_findClusterMates_write_log;
#endif
template<typename TCdrGlobalData>
void findClusterMates(TCloneStore::value_type const & queryElement, 
        TCloneStore::const_iterator const & beginRefElements, 
        TCloneStore::const_iterator const & endRefElements,
        TCdrGlobalData & global,
        TLQPositionStore & lqPositionStore, 
        CloneToPairsMap & clusterPairsByMinor,
        ClusterStats & clusterStats,
        ProgressBar & progBar) 
{
    CdrOptions const & options = global.options;

    unsigned nPwAlignments = 0;
    unsigned long long counter = 0;
    unsigned localTrueNegatives = 0, localFalseNegatives = 0;

    std::stringstream clusterEvalLogSS;

    for (TCloneStore::const_iterator tar = beginRefElements; tar!=endRefElements; ++tar) {
        ++counter;

        TCloneStore::value_type const & targetElement = *tar;

        double countRatio = static_cast<double>(queryElement.second.count) / static_cast<double>(targetElement.second.count);
        if (countRatio > 1)
            countRatio = 1.0 / countRatio;

        ClusterFailure hasFailed = NO_FAILURE;

        // Check if the two CDR3 sequences have the same length
        unsigned lenDiff = std::abs(static_cast<int>(length(queryElement.first.cdrSeq)) - static_cast<int>(length(targetElement.first.cdrSeq)));
        if (lenDiff) {
            hasFailed = CDR3LEN_MISMATCH;
        }

        // Check if the two sets of V alleles and the two sets of J alleles
        // intersect
        else if ( !doIntersectSTL(queryElement.first.VIds, targetElement.first.VIds) || !doIntersectSTL(queryElement.first.JIds, targetElement.first.JIds)) {
            hasFailed = VJ_MISMATCH;
        } 
        
        // Check if the ratio below or equal to the user specified threshold
        else if (countRatio > global.options.maxClusterRatio) {
            hasFailed = RATIO_EXCEEDED;
        }

        else {
            unsigned maxQClusterErrors = options.qualClustering ? options.maxQClusterErrors : 0;
            unsigned maxSClusterErrors = options.simpleClustering ? options.maxSClusterErrors : 0;
            unsigned maxErrors = maxQClusterErrors + maxSClusterErrors;

            ClusterPair cp(queryElement, targetElement, maxErrors);
            ++nPwAlignments;

            // Only if the alignment of the two clusters contains at most maxQClusterErrors errors
            // the cluster pair is stored
            if (!cp.isValid()) {
                hasFailed = TOO_MANY_ERRORS;
            }

            unsigned nErrors = length(cp.getMinorErrPositions());
            // Theoretically we don't have to check the quality scores of the
            // error positions if #errors is below the simple clustering
            // threshold, but we do it for the stats
            unsigned lqCorrelated = 0;
            if (options.qualClustering && (!hasFailed)) {
                // Check if all error positions correspond to a low quality score in the minor
                // clonotype
                Clone<Dna5> c = cp.getMinorClone();
                TCloneStore::value_type const & minorElement = (c==queryElement.first) ? queryElement : targetElement;
                std::set<unsigned> minorLQPositions = getLQPositions(lqPositionStore, minorElement.second, options);
                for (Iterator<String<unsigned> const ,Rooted>::Type minErr=begin(cp.getMinorErrPositions()); !atEnd(minErr); goNext(minErr)) {
                    if (minorLQPositions.find(*minErr) != minorLQPositions.end()) {
                        ++lqCorrelated;
                    }
                }
                // If we have more LQ correlated errors than allowed, we drop
                // them so that they have to be covered by the simple
                // clustering errors
                lqCorrelated = std::min(maxQClusterErrors, lqCorrelated);
                SEQAN_CHECK(nErrors >= lqCorrelated, "PLEASE REPORT THIS ERROR");
                if ((nErrors - lqCorrelated) > maxSClusterErrors) {
                    hasFailed = TOO_MANY_ERRORS;
                }
            }

            // If all checks passed, store the clusterpair. Not thread-safe, has
            // to be protected.
            if (!hasFailed) {
                if (lqCorrelated == 0)
                    cp.reason = "SCLUST";
                else if (nErrors == lqCorrelated)
                    cp.reason = "QCLUST";
                else
                    cp.reason = "QCLUST+SCLUST";
#ifdef __WITHCDR3THREADS__
                std::lock_guard<std::mutex> lock(MUTEX_findClusterMates_store_result);
#endif
                clusterPairsByMinor[cp.getMinorClone()].insert(cp);
                continue; // <= Code execution ends here for POSITIVES
            } 
        }

    }

#ifdef __WITHCDR3THREADS__
    std::lock_guard<std::mutex> lock(MUTEX_findClusterMates_write_log);
#endif

    if (global.outFiles.clusterEvalLog.wasSpecified())
        global.outFiles.clusterEvalLog.getStream() << clusterEvalLogSS.str();

    clusterStats.alignmentCount += nPwAlignments;
    clusterStats.checkedPairs += counter;
    clusterStats.falseNegatives += localFalseNegatives;
    clusterStats.trueNegatives += localTrueNegatives;
    progBar.updateAndPrint(counter);
}

inline std::string offOnBool(bool const b) {
    if (b)
        return "on";
    else
        return "off";
}

template<typename TCdrGlobalData>
inline void runClonotypeClustering(TCdrGlobalData & global, Dna5CloneStore & cloneStore) 
{
    // Check if any clustering mode was enabled
    if ((!global.options.simpleClustering) && (!global.options.qualClustering))
        return;

    typedef std::multiset<const Clone<Dna5>*, ClonesBySize > TClonePointerMultiSet;

    std::cerr << "===== Posterior clustering of clonotypes" << std::endl;
    std::cerr << "  |-- Quality clustering: " << offOnBool(global.options.qualClustering) << "\n";
    std::cerr << "  |-- Simple clustering: " << offOnBool(global.options.simpleClustering) << "\n";

    // ============================================================================
    // Compute all pairwise alignments and filter the potential major clonotypes
    // ============================================================================

    unsigned long nPairs = cloneStore.size();
    nPairs = (nPairs * nPairs - nPairs) / 2;

    std::cerr << "  |-- Number of clone store elements: " << cloneStore.size() << std::endl;
    std::cerr << "  |-- Number of alignment pairs: " << nPairs<< std::endl;

    CloneToPairsMap clusterPairsByMinor(cloneStore);

    std::cerr << "  |-- Computing alignments" << std::endl;

    TLQPositionStore lqPositionStore;

    ProgressBar progBar(std::cerr, nPairs, 100, "      ");
    progBar.print_progress();

    ClusterStats clusterStats;

    // -!- ==========================================================================
    // -!- Multithreading dependent code - if compiled with multi-threading support,
    // -!- a ThreadPool is created and handles one task for every for-loop
    // -!- iteration.  Otherwise, the for-loop does not enqueue tasks in the thread
    // -!- pool but actually executes the function.
    // -!- ==========================================================================
    std::clock_t clockBeforeClustAlign = std::clock();
    {

#ifdef __WITHCDR3THREADS__
        ThreadPool threadPool(global.options.jobs);
#endif

        unsigned __cnt = 0;
        for (Dna5CloneStore::const_iterator cluster1 = cloneStore.begin(); cluster1 != cloneStore.end(); ++cluster1, ++__cnt) {
            Dna5CloneStore::const_iterator cluster2 = cluster1;
            ++cluster2;
#ifdef __WITHCDR3THREADS__
            threadPool.enqueue<void>([cluster1, cluster2, &cloneStore, &global, &lqPositionStore, &clusterPairsByMinor, &clusterStats, &progBar]()
                    {
#endif
                    findClusterMates(*cluster1, cluster2, cloneStore.end(), global, lqPositionStore, clusterPairsByMinor, clusterStats, progBar);
#ifdef __WITHCDR3THREADS__
                    }
                    );
#endif

        }
    } // Destructs ThreadPool, joins all threads

    progBar.clear();

    std::cerr << "  |-- Required cpu time: " << formatSeconds(double(std::clock() - clockBeforeClustAlign) / CLOCKS_PER_SEC) << std::endl;;
    std::cerr << "  |-- computed " << clusterStats.alignmentCount << " pairwise alignments." << std::endl;
//    std::cerr << "  |~~ [DEBUG] Number of checked pairs: " << clusterStats.checkedPairs << std::endl;

    // ==========================================================================
    // Store the minors by size
    // ==========================================================================

    ClonesBySize cbs(cloneStore);

    TClonePointerMultiSet minorsBySize(cbs);

    for (CloneToPairsMap::TStore::const_iterator it = clusterPairsByMinor.getMap().begin(); it!=clusterPairsByMinor.getMap().end(); ++it)
        minorsBySize.insert(&(it->first));

    // ==========================================================================
    // DEBUG
    // ==========================================================================

    uint64_t foundPairs = 0;
    for (CloneToPairsMap::TStore::const_iterator it = clusterPairsByMinor.getMap().begin(); it!=clusterPairsByMinor.getMap().end(); ++it)
        foundPairs += it->second.size();

    // ============================================================================
    // Output some info about the potential clustering impact
    // ============================================================================

    unsigned nAnalyzedReads = 0;
    for (Dna5CloneStore::const_iterator cse = cloneStore.begin(); cse!=cloneStore.end(); ++cse)
        nAnalyzedReads += cse->second.count;
    unsigned nReadsAffected = 0;
    unsigned nGoodPairs = 0;
    for (CloneToPairsMap::TStore::const_iterator cp = clusterPairsByMinor.getMap().begin(); cp!=clusterPairsByMinor.getMap().end(); ++cp) {
        nGoodPairs += length(cp->second);
        nReadsAffected += cloneStore[cp->first].count;
    }

    std::cerr << "  |-- " << nGoodPairs << " cluster pairs with at most " << global.options.maxQClusterErrors << " errors were found." << std::endl;
    std::cerr << "  |-- " << std::setprecision(2) << std::fixed << nReadsAffected * 100.0 / nAnalyzedReads 
        << " % of the " << nAnalyzedReads 
        << " successfully analyzed reads are similar to another\n      cluster within the specified error margin." << std::endl;;
    std::cerr << "  |-- " << std::setprecision(2) << std::fixed << 100.0 * clusterPairsByMinor.getMap().size() / cloneStore.size() << " % of the " 
        << cloneStore.size() << " identified clonotypes are potentially affected by clustering" << std::endl;

    // ============================================================================
    // Iterate through all minors from the smallest to the larges and check for
    // each if it can be corrected
    // ============================================================================

    ClusterLog clusterLog;

    for (TClonePointerMultiSet::const_iterator it = minorsBySize.begin(); it!=minorsBySize.end(); ++it) {
        CloneToPairsMap::TStore::mapped_type & majors = clusterPairsByMinor[**it];
        Clone<Dna5> modClone = **it;
        Clone<Dna5> const & minorClone = **it;

        global.outFiles.clusterCLog <<minorClone;

        ClusterResult const & minorResult = cloneStore[minorClone];

        global.outFiles.clusterCLog << "\t";
        unsigned maxNErrors = 0;
        unsigned i=1;
        for (CloneToPairsMap::TStore::mapped_type::const_iterator m = majors.begin(); m!=majors.end(); ++m)  {
            if (length(m->getMinorErrPositions()) > maxNErrors)
                maxNErrors = length(m->getMinorErrPositions());
            global.outFiles.clusterCLog << i++ << ":" << m->getMajorClone() << "(" << m->reason << ");";
        }

        String<ClusterPair> targetPairs;
        for (CloneToPairsMap::TStore::mapped_type::const_iterator m = majors.begin(); m!=majors.end(); ++m) 
            if (length(m->getMinorErrPositions()) == maxNErrors)
                appendValue(targetPairs, *m);

        String<unsigned> redistributionCounts;
        String<double> redistributionRatios;
        computeRedistCounts(redistributionCounts, redistributionRatios, minorClone, targetPairs, cloneStore);
        global.outFiles.clusterCLog << "\t";
        for (Iterator<String<unsigned>,Rooted>::Type it = begin(redistributionCounts); !atEnd(it); goNext(it))
            global.outFiles.clusterCLog << *it << ";";
        global.outFiles.clusterCLog << std::endl;
//        ::operator<<(global.outFiles.clusterCLog << "REDIST\t" ,  redistributionCounts) << std::endl;

        // ============================================================================
        // Iterate trough all the major CTs and "upgrade" them according to the 
        // computed redistribution counts
        // ============================================================================

        String<Clone<Dna5> > majorClones;
        for (Iterator<String<ClusterPair>, Rooted>::Type cpIt = begin(targetPairs); !atEnd(cpIt); goNext(cpIt)) {
            ClusterResult partialMinorResult = minorResult;
            partialMinorResult.count = redistributionCounts[position(cpIt)];
            Dna5CloneStore::iterator cseMajor = cloneStore.find(cpIt->getMajorClone());
            if (cseMajor == cloneStore.end()) {
                std::cerr << "[ERR] Error in runClonotypeClustering(); cannot find major clone" << std::endl;
                exit(1);
            }

            // Store the major if we are logging
            if (global.outFiles.clusterEvalLog.good())
                appendValue(majorClones, cseMajor->first);

            // Do the actual clustering
            mergeWithClusterResult(cseMajor->second, partialMinorResult, cpIt->getMinorErrPositions());
        }

        // ============================================================================
        // If enabled, log this cluster event
        // ============================================================================
        if (global.outFiles.clusterEvalLog.good())
            clusterLog.logClusterEvent(minorClone, majorClones, cloneStore, redistributionRatios);

        // ============================================================================
        // After the read count of the minor CT has been redistributed entirely, delete
        // the clone
        // ============================================================================

        Dna5CloneStore::iterator testIt = cloneStore.find(minorClone);
        if (testIt == cloneStore.end())
            exitWithError("ERROR 9999 - didnt find clone for delection", 9999);

        cloneStore.erase(testIt);

    }

}

        

/*
 * Translates clonotype counts from nucleotide based clonotypes to aminoacid based
 * clonotype counts
 */
template <typename TNucAcid>
void _translateClones(std::map<Clone<AminoAcid>, ClusterResult>& aaClones, std::map<Clone<TNucAcid>, ClusterResult> const & nucClones) {
    aaClones.clear();
    
    for (typename std::map<Clone<TNucAcid>, ClusterResult>::const_iterator nucClone = nucClones.begin(); nucClone!=nucClones.end(); ++nucClone) {
        Clone<AminoAcid> aaClone;
        aaClone.VIds = nucClone->first.VIds;
        aaClone.JIds = nucClone->first.JIds;
        translate(aaClone.cdrSeq, nucClone->first.cdrSeq);
        aaClones[aaClone].count += nucClone->second.count;
        //TODO merge average qualities
    }
}

/*
 * Writes clonotype counts to a stream
 */
template <typename TAlphabet,typename TCdrGlobalData>
void _writeClonotypeCounts(std::ostream& streamCounts, std::map<Clone<TAlphabet>, ClusterResult>const & cloneStore, TCdrGlobalData & global, bool translateToAA = false)  {

    typedef    std::map<CharString, unsigned> TCounter;

    bool mergeAllels = true;
    TCounter fingerPrintCounter;

    // ============================================================================
    // Build the fingerprint strings for all clones and increase the counter
    // ============================================================================

    for (typename std::map<Clone<TAlphabet>,ClusterResult>::const_iterator storeElem = cloneStore.begin(); storeElem!=cloneStore.end(); ++storeElem) {
        const Clone<TAlphabet>& clone = storeElem->first;

        CharString fingerPrint;
        _appendGeneList(fingerPrint, clone.VIds, global.references.leftMeta, mergeAllels);
        appendValue(fingerPrint, ':');
        append(fingerPrint, clone.cdrSeq);
        appendValue(fingerPrint, ':');
        _appendGeneList(fingerPrint, clone.JIds, global.references.rightMeta, mergeAllels);
        if (translateToAA) {
            String<AminoAcid> translatedCDR;
            translate(translatedCDR, clone.cdrSeq);
            appendValue(fingerPrint, '\t');
            append(fingerPrint, translatedCDR);
        }
        fingerPrintCounter[fingerPrint] += storeElem->second.count;
    }

    // ============================================================================
    // Write the counts to the output streamCounts
    // ============================================================================
    for (TCounter::const_iterator fpCount = fingerPrintCounter.begin(); fpCount!=fingerPrintCounter.end(); ++fpCount)
        streamCounts << fpCount->first << '\t' << fpCount->second << std::endl;

}

void closeOfStream(std::ofstream *s) {
    if (s!=NULL) {
        s->close();
        delete s;
    }
}

void initOutFileStream(CharString const & path, std::ofstream*& target) throw(std::string) {
    if (path!="") {
        target = new std::ofstream(toCString(path));
        if (!target->good())
            throw std::string("Cannot open " + std::string(toCString(path)));
    } else {
        target = NULL;
    }
}

template<typename TRow>
void writeAlign(std::ostream& out, TRow& row0, TRow& row1) {
    Align<typename Source<TRow>::Type> align;
    appendValue(align.data_rows, row0);
    appendValue(align.data_rows, row1);
    out << align;
}

template<typename TRow>
inline unsigned viewEndPosition(TRow const & row) {
    return toViewPosition(row, endPosition(row));
}

struct ErrorStats {
    unsigned matchLen, outerMatchLen;
    double errRate, outerErrRate;
    inline std::string toString() const;
    String<int> outerErrPositions;

    std::string toString() {
        std::stringstream ss;
        ss << "matchLen = " << matchLen << "; outerMatchLen = " << outerMatchLen << "; errRate = " << errRate
            << "; outerErrRate = " << outerErrRate << "; outerErrPositions = {";
        for (auto & x : outerErrPositions)
            ss << x << " ";
        ss << "}";
        return ss.str();
    }
};

inline std::string ErrorStats::toString() const {
    std::stringstream ss;
    ss << "{matchLen=" << matchLen << "; outerMatchLen=" << outerMatchLen << "; errRate=" << errRate << "; outerErrRate=" << outerErrRate << "}";
    return ss.str();
}

template<typename TRow, typename TPos>
inline bool noGap(TRow const & row1, TRow const & row2, TPos pos) {
    bool r1 = isGap(row1, pos);
    bool r2 = isGap(row2, pos);
    SEQAN_CHECK(!(r1 && r2), "ERROR 1011 - double gap position in pairwise alignment [noGap()]. Please report this error!");
    return !r1 && !r2;
}

template<typename TAlign, typename TPos>
inline bool noGap(TAlign const & align, TPos pos) {
    bool r1 = isGap(row(align, 0), pos);
    bool r2 = isGap(row(align, 1), pos);
    if (r1 && r2) {
        std::cerr << align << std::endl;
        std::cerr << "CRITICAL POSITION: " << pos << std::endl;
    }

    SEQAN_CHECK(!(r1 && r2), "ERROR 1011 - double gap position in pairwise alignment [noGap()]. Please report this error!");
    return !r1 && !r2;
}

/**
 * Counts the number of matches outside the CDR3 region plus the number of matches
 * at most N bases into the CDR3 region, where N is a user specified parameter
 */
template<typename TSeq>
ErrorStats computeErrorStats (
        Align<TSeq> const & align,                  // An alignment read vs. reference segment
        SegmentMeta const & meta,            // The meta information for the aligned reference
        LeftOverlap const &)                 // Tag indicating the orientation of the overlap <RightOverlap|LeftOverlap>
{

    typedef typename Row<Align<TSeq> const>::Type    TRow;
    typedef typename Position<TRow>::Type       TPos;

    TRow& readRow       = row(align, 0);
    TRow& segmentRow    = row(align, 1);

    // -----------xxxxxxxxxxxxxCysxxxxxxxxxxxxxxxxxxxxxxxxx     [read]
    // xxxxxxxxxxxxxxxxxxxxxxxxCysxxxxxxxxxxxxx------------     [V-segment]
    // 0          begin           end         maxEnd

    TPos maxEnd = std::min(viewEndPosition(readRow), viewEndPosition(segmentRow));

    TPos begin = std::max(toViewPosition(readRow, 0), toViewPosition(segmentRow, 0));
    // TODO Check if motif position is covered by the alignment. In principle,
    // this is guaranteed by the positioning of the segment core fragment
    TPos x = meta.motifPos - beginPosition(source(segmentRow));
    TPos y = toViewPosition(segmentRow, x);
    TPos end = std::min(y + 3, maxEnd);

    ErrorStats es;
    unsigned nMatches = 0, nTotalMatches = 0;

    TPos pos = begin;
    for (; pos < end; pos++) {
        if (noGap(readRow, segmentRow, pos) && readRow[pos] == segmentRow[pos])
            ++nMatches;
        else
            appendValue(es.outerErrPositions, -(toSourcePosition(segmentRow, end) - toSourcePosition(segmentRow, pos))+1);
    }
    for (; pos < maxEnd; pos++)
        if (noGap(readRow, segmentRow, pos) && readRow[pos] == segmentRow[pos])
            ++nTotalMatches;
    nTotalMatches += nMatches;

    es.outerErrRate = 1.0 - (double(nMatches) / (end-begin));
    es.errRate = 1.0 - (double(nTotalMatches) / (maxEnd-begin));
    es.matchLen  = maxEnd - begin;
    es.outerMatchLen = end - begin;

    return es;
}

template<typename TSeq>
ErrorStats computeErrorStats (
        Align<TSeq> const & align,       // An alignment read vs. reference segment
        SegmentMeta const & meta,        // The meta information for the aligned reference
        RightOverlap const &)            // Tag indicating the orientation of the overlap <RightOverlap|LeftOverlap>
{
    typedef typename Row<Align<TSeq> const>::Type    TRow;
    typedef typename Position<TRow>::Type       TPos;

    TRow& readRow       = row(align, 0);
    TRow& segmentRow    = row(align, 1);
    
    TPos maxEnd = std::min(viewEndPosition(readRow), viewEndPosition(segmentRow));
    TPos minBegin = std::max(toViewPosition(readRow, 0), toViewPosition(segmentRow, 0));
    // TODO Check if motif position is covered by the alignment. In principle,
    // this is guaranteed by the positioning of the segment core fragment
    TPos motifViewPosition = toViewPosition(segmentRow, meta.motifPos - beginPosition(source(segmentRow)));

    TPos begin = std::max(minBegin, motifViewPosition);
    //                    ^- is this necessary?

    // xxxxxxxxxxxxxxxxxxxxxxxxPhexxxxxxxxxxxxx------------     [read]
    // -----------xxxxxxxxxxxxxPhexxxxxxxxxxxxxxxxxxxxxxxxx     [J-segment]
    // 0          minBegin     begin           maxEnd

    unsigned nMatches = 0, nTotalMatches = 0;
    ErrorStats es;
    TPos pos = minBegin;
    for (; pos <= begin; pos++) 
        if (readRow[pos] == segmentRow[pos])
            ++nTotalMatches;
    for (; pos < maxEnd; pos++)  {
        if (readRow[pos] == segmentRow[pos])
            ++nMatches;
        else
            appendValue(es.outerErrPositions, toSourcePosition(segmentRow, pos) - toSourcePosition(segmentRow, begin));
    }
    nTotalMatches += nMatches;

    es.matchLen = maxEnd - (minBegin);
    es.outerMatchLen = maxEnd - begin;
    es.outerErrRate = es.outerMatchLen == 0 ? 1.0 : 1.0 - double(nMatches) / es.outerMatchLen;
    es.errRate = es.matchLen == 0 ? 1.0 : 1.0 - double(nTotalMatches) / es.matchLen;

    return es;
}

template<typename TSeq, typename TOverlapDirection>
bool getErrorPositions(
        String<int> & errorPositions,                  // [OUT] The error positions
        unsigned & outerMatchLen,                      // [OUT] The length of the outer alignment
        String<SegmentMatch<TSeq> > const & matches,   // [IN]  The best left segment matches
        String<SegmentMeta> const & meta,              // [IN]  The meta information for the aligned reference
        TOverlapDirection const &)                     // [TAG] Indicating the orientation of the overlap <RightOverlap|LeftOverlap>
{
    clear(errorPositions);
    for (unsigned i = 0; i<length(matches); ++i) {
        ErrorStats es = computeErrorStats(matches[i].align, meta[matches[i].db], TOverlapDirection());
        if (i==0) {
            errorPositions = es.outerErrPositions;
            outerMatchLen  = es.outerMatchLen;
        } else if (es.outerErrPositions != errorPositions) {
            clear(errorPositions);
            outerMatchLen  = -1u;
            return false;
        }
    }
    return true;
}

/*
 * Corrects the alignment score taking into account the position of the
 * left/right motif.  Positions within the CDR3 do not contribute negatively to
 * the score, as otherwise longer segments would loose more score than shorter
 * segments, give RAG cut rather early. Positive contributions are allowed if
 * uninterrupted from the beginning of the CDR3.
 */
template<typename TSeq>
inline void _refineAlignmentScore(
        int & scoreVal,                      // The alignment score previously computed
        Align<TSeq>& align,                  // An alignment read vs. reference segment
        SegmentMeta const & meta,            // The meta information for the aligned reference
        Score<int> const & scoringScheme,    // The scoring scheme for the candidate verification alignment
        LeftOverlap const &)                 // Tag indicating the orientation of the overlap <RightOverlap|LeftOverlap>
{
    typedef typename Row<Align<TSeq> >::Type     TRow;

    TRow& readRow       = row(align, 0);
    TRow& segmentRow    = row(align, 1);

    // Checked, still works with new align module
    unsigned begin = toViewPosition(segmentRow, meta.motifPos) + 3;
    unsigned end = std::min(viewEndPosition(readRow), viewEndPosition(segmentRow));

    bool tail = true;
    for (unsigned pos = begin;pos < end; ++pos) {
        if (tail && readRow[pos] == segmentRow[pos])
            continue;
        tail = false;
        scoreVal -= score(scoringScheme, readRow[pos], segmentRow[pos]);
    }

}

//TODO alignment score refinement
//
// =1= Truncate all V/J sequences based on the user specified v- or j-boundary.
// =2= Make sure that all V/J segments are evaluated according to the same base
//     length, i.e. that the tails have the same lengths. Check if there are any
//     relevant cases in the reference files.

template<typename TSeq>
inline void _refineAlignmentScore(
        int & scoreVal,                      // The alignment score previously computed
        Align<TSeq>& align,                  // An alignment read vs. reference segment
        SegmentMeta const & meta,            // The meta information for the aligned reference
        Score<int> const & scoringScheme,    // The scoring scheme for the candidate verification alignment
        RightOverlap const &)                // Tag indicating the orientation of the overlap <RightOverlap|LeftOverlap>
{
    typedef typename Row<Align<TSeq> >::Type     TRow;

    TRow& readRow       = row(align, 0);
    TRow& segmentRow    = row(align, 1);
    int pos = toViewPosition(segmentRow, meta.motifPos) - 1;
    int end = std::max(beginPosition(readRow), beginPosition(segmentRow));
    //int end = std::max(toViewPosition(readRow, 0), toViewPosition(segmentRow, 0));

    bool tail = true;
    for (;pos >= end; --pos) {
        if (tail && readRow[pos] == segmentRow[pos])
            continue;
        tail = false;
        scoreVal -= score(scoringScheme, readRow[pos], segmentRow[pos]);
    }

}

template<typename TSeq, typename TOriSpec>
inline void _computeAlignMeta(SegmentMatch<TSeq> & match, SegmentMeta const & meta, TOriSpec const &) {

    typedef typename Row<Align<TSeq> >::Type     TRow;
    
    TRow& readRow       = row(match.align, 0);
    TRow& segmentRow    = row(match.align, 1);

    match.matchPhreds.clear();
    match.mismatchPhreds.clear();

    // ============================================================================
    // Compute the boundaries of the non-cdr3 region of interest
    // ============================================================================

    int motifPos = toViewPosition(segmentRow, meta.motifPos);
    int start,end;

    if (IsSameType<TOriSpec, RightOverlap>::VALUE) {
        start = motifPos;
        end   = std::min(viewEndPosition(readRow), viewEndPosition(segmentRow));
    } else if (IsSameType<TOriSpec, LeftOverlap>::VALUE) {
        start = std::max(toViewPosition(readRow, 0), toViewPosition(segmentRow, 0));
        end   = motifPos + 3;
    }

    // ============================================================================
    // Count the errors within the non-cdr3 alignment
    // ============================================================================

    unsigned errors = 0;
    match.outerErrPos.clear();
    match.motifPos = toSourcePosition(readRow, motifPos);

    for (int pos = start; pos<end; ++pos) {
        if (isGap(readRow,pos) || isGap(segmentRow, pos) || segmentRow[pos] != readRow[pos]) {
            ++errors;
            // Store the read error positions
            match.outerErrPos.push_back(toSourcePosition(readRow, pos));
            match.mismatchPhreds[getQualityValue(convert<Dna5Q>(readRow[pos]))]++;
        } else {
            match.matchPhreds[getQualityValue(convert<Dna5Q>(readRow[pos]))]++;
        }
    }

    // ============================================================================
    // Compute the error rate and alignment length
    // ============================================================================
    
    match.outerLen      = end - start;
    match.outerErrRate  = errors / ( 1.0 * match.outerLen );

    // ============================================================================
    // Write the mismatch counts for each score
    // ============================================================================

    PHREDSTAT << source(readRow) << "," << match.outerErrRate;

    if (IsSameType<TOriSpec, RightOverlap>::VALUE) {
        PHREDSTAT << ",R";
    } else if (IsSameType<TOriSpec, LeftOverlap>::VALUE) {
        PHREDSTAT << ",L";
    } else {
        std::cerr << "ERR" << std::endl;
        SEQAN_ASSERT_MSG(false,"INVALID TAG");
    }

    PHREDSTAT << ",MISMATCH";

    for (unsigned score = 0; score <= 60; ++score)
        PHREDSTAT << "," << match.mismatchPhreds[score];

    PHREDSTAT << std::endl;

    // ============================================================================
    // Write the match counts for each score
    // ============================================================================

    PHREDSTAT << source(readRow) << "," << match.outerErrRate;

    if (IsSameType<TOriSpec, RightOverlap>::VALUE) {
        PHREDSTAT << ",R";
    } else if (IsSameType<TOriSpec, LeftOverlap>::VALUE) {
        PHREDSTAT << ",L";
    } else {
        std::cerr << "ERR" << std::endl;
        SEQAN_ASSERT_MSG(false,"INVALID TAG");
    }

    PHREDSTAT << ",MATCH";

    for (unsigned score = 0; score <= 60; ++score)
        PHREDSTAT << "," << match.matchPhreds[score];

    PHREDSTAT << std::endl;

}

template <typename TSequence>
unsigned getAlignmentErrors(Align<TSequence> const & align, unsigned segBeginPos, unsigned segEndPos) {
    typedef Align<TSequence> const          TAlign;
    typedef typename Row<TAlign>::Type      TRow;
    typedef typename Source<TRow>::Type     TSource;

    TRow& readRow = row(align,0);
    TRow& segRow  = row(align,1);

    TSource & segSource = source(segRow);

    unsigned viewEndPos = std::min(viewEndPosition(readRow), viewEndPosition(segRow));

    // Adjust the passed begin and end positions in case the source of the
    // alignment rows is a 'Segment'
    segBeginPos -= beginPosition(segSource);
    segEndPos -= beginPosition(segSource);
    
    // If an end position behind the source sequence was specified, fix that
    segEndPos = std::min(segEndPos, static_cast<unsigned>(length(source(segRow))));

    if (segEndPos < segBeginPos) {
            std::cerr << "[ERR] end < begin occurred in getAlignmentErrors(). Please report this error." << std::endl;
            exit(1);
    }

    unsigned errors = segEndPos - segBeginPos;
    for (unsigned pos=toViewPosition(segRow, segBeginPos); pos < viewEndPos && pos<=toViewPosition(segRow, segEndPos-1) ; ++pos) {
    // This breaks when checking for segEndPos, if segEndPos is tight.                                     ^^^^^^^^^^^
    // Therefore, check for <=segEndPos-1
        if (!isGap(readRow, pos) && readRow[pos] == segRow[pos])
            --errors;
    }

    return errors;
}

template <typename TSequence>
inline unsigned getReadMotifPos(SegmentMatch<TSequence> const & match, String<SegmentMeta> const & meta) {
    // Type definitions
    typedef typename SegmentMatch<TSequence>::TAlign     TAlign;
    typedef typename Row<TAlign>::Type                              TRow;
    typedef typename Position<TSequence>::Type                      TPos;

    TRow const & readRow = row(match.align, 0);
    TRow const & segRow = row(match.align, 1);

    // This takes into account that source of the alignment rows can
    // potentially be 'Segment's. beginPosition() and endPosition() return 0
    // for 'String' types and therefore don't influence the result if this is
    // not the case. 
    TPos segInfixMotifPos = meta[match.db].motifPos - beginPosition(source(segRow));
    return toSourcePosition(readRow, toViewPosition(segRow, segInfixMotifPos)) + beginPosition(source(readRow));
}

template <typename TInfixHost, typename TSequence, typename TGlobal>
RejectReason findCDR3Region(
        Segment<TInfixHost, InfixSegment> & infix,                // The Infix<> object to hold the output substring
        String<SegmentMatch<TSequence> > const & leftMatches,   // The best left segment matches
        String<SegmentMatch<TSequence> >  const & rightMatches,  // The best right segment matches
        TGlobal & global)                                // The CdrOptions
{
    typedef String<SegmentMatch<TSequence> > const    TSegMatches;
    typedef typename Iterator<TSegMatches, Rooted>::Type        TSegmIterator;

    if (length(leftMatches)<1 || length(rightMatches)<1 ) {
        std::cerr << "\n\nERROR: Attempted to identify CDR3 region without V or J matches. Please report this error!" << std::endl;
        exit(1);
    }

    unsigned leftMotifReadPos = getReadMotifPos(leftMatches[0], global.references.leftMeta);
    unsigned rightMotifReadPos = getReadMotifPos(rightMatches[0], global.references.rightMeta);

    // ============================================================================
    // If we have more than one match at any of the sides, make sure that all
    // candidates have their motif position aligned to the same location within
    // the read. Otherwise the CDR3 begin / end position is ambiguous and the read
    // is rejected.
    // ============================================================================

    TSegmIterator leftIt = begin(leftMatches);
    for (goNext(leftIt); !atEnd(leftIt); goNext(leftIt)) {
        if (leftMotifReadPos != getReadMotifPos(*leftIt, global.references.leftMeta)) 
            return MOTIF_AMBIGUOUS;
    }

    TSegmIterator rightIt = begin(rightMatches);
    for (goNext(rightIt); !atEnd(rightIt); goNext(rightIt))
        if (rightMotifReadPos != getReadMotifPos(*rightIt, global.references.rightMeta)) 
            return MOTIF_AMBIGUOUS;

    setBeginPosition(infix, leftMotifReadPos);
    setEndPosition(infix, rightMotifReadPos+3);

    return NONE;
}

const AnalysisResult rejectRead(String<Dna5> const & readSeq, RejectReason const & reason, TCacheStore* cacheStore) {
    AnalysisResult as(reason);
    CacheEntry ce(as);
#ifdef __WITHCDR3THREADS__
    std::unique_lock<std::mutex> lock(CACHE_STORE_MUTEX);
#endif
    //std::cerr << "REJECT FOR " << reason << std::endl;
    if (cacheStore != NULL) {
        std::pair<TCacheStore::iterator,bool> prev = cacheStore->insert( std::make_pair(readSeq, ce) );
        if (!prev.second && prev.first->second != ce) {
            std::cerr << "\n\nERROR: Cache consistency error [reject]" << std::endl;
            exit(1);
        }
    }
    return as;
}

template<typename TContainer>
inline void getQualityString(String<uint64_t>& targetString, TContainer const & dnaString) {
    resize(targetString, length(dnaString));
    for (typename Iterator<TContainer const, Rooted>::Type ch = begin(dnaString); !atEnd(ch); goNext(ch)) 
        targetString[position(ch)] = getQualityValue(*ch);
}

inline void countNewClone(TCloneStore & clusterStore, Clone<Dna5> const & clone, String<uint64_t> const & qualities) {
    ClusterResult& result = clusterStore[clone];
    result.count++;
    if (result.count == 1) { // First clone of this kind
        result.sumOfQScores = qualities;
    } else {
        if (length(qualities) != length(result.sumOfQScores)) {
            std::cerr << "ERROR: Quality vector lengths don't match. Please report this error!" << std::endl;
            exit(1);
        }
        for (Iterator<String<uint64_t>, Rooted>::Type qIt = begin(result.sumOfQScores); !atEnd(qIt); goNext(qIt))  
            (*qIt) += qualities[position(qIt)];

        //TODO don't divide for each newly encountered read but only once at the end

    }
}

/**
 * Truncate the sequences in a FastqRecord - truncating always happens from the
 * *end* of each read since it is performed before any reverse-complement
 * is built!
 *
 * The single end implementation ignores the vReadCrop parameter.
 */
inline void truncateRead(FastqRecord<SingleEnd> & fqRecord, unsigned len, unsigned vReadCrop) {
    ignoreUnusedVariableWarning(vReadCrop);
    if (len > 0 && length(fqRecord.seq) > len)
        resize(fqRecord.seq, len);
}

inline void truncateRead(FastqRecord<PairedEnd> & fqRecord, unsigned len, unsigned vReadCrop) {
    if (len > 0 && length(fqRecord.fwSeq) > len)
        resize(fqRecord.fwSeq, len);
    if (vReadCrop > 0) {
        if (length(fqRecord.fwSeq) > vReadCrop)
            erase(fqRecord.fwSeq, 0, vReadCrop);
        else
            clear(fqRecord.fwSeq);
    }
    if (len > 0 && length(fqRecord.revSeq) > len)
        resize(fqRecord.revSeq, len);
}

///*
// * Returns true if the average quality of the reverse read (paired end) or the read
// * (single end) is above or equal the user specified threshold
// */
//inline bool checkQuality(FastqRecord<SingleEnd> & fqRecord, CdrOptions const & options) 
//{
//    double qualAvg = 0;
//    for (unsigned pos=0; pos<length(fqRecord.seq); ++pos)
//        qualAvg += getQualityValue(fqRecord.seq[pos]);
//    qualAvg /= length(fqRecord.seq);
//
//    return (qualAvg >= options.qmin);
//}
//
//inline bool checkQuality(FastqRecord<PairedEnd> & fqRecord, CdrOptions const & options) 
//{
//    // If the reverse (VDJ) read doesn't pass the QC, reject the record
//    double qualAvg = 0;
//    for (unsigned pos=0; pos<length(fqRecord.revSeq); ++pos)
//        qualAvg += getQualityValue(fqRecord.revSeq[pos]);
//    qualAvg /= length(fqRecord.revSeq);
//
//    if (qualAvg < options.qmin)
//        return false;
//
//    // If the forward (V) read doesn't pass the QC, clear the sequence so
//    // that a V-identification can still be attempted on the reverse (VDJ) read.
//    qualAvg = 0;
//    for (unsigned pos=0; pos<length(fqRecord.fwSeq); ++pos)
//        qualAvg += getQualityValue(fqRecord.fwSeq[pos]);
//    qualAvg /= length(fqRecord.fwSeq);
//
//    if (qualAvg < options.qmin)
//        clear(fqRecord.fwSeq);
//
//    return true;
//}

inline bool inStreamsAtEnd(SeqInputStreams<PairedEnd> const & inStreams) 
{
    bool s1End = atEnd(inStreams.fwStream);
    bool s2End = atEnd(inStreams.revStream);

    if (s1End && s2End)
        return true;
    if (s1End || s2End) {
        std::cerr << "\nWARNING - Reached end of one input file before reaching the end of the other!\n";
        return true;
    }
    return false;
}

inline bool inStreamsAtEnd(SeqInputStreams<SingleEnd> const & inStreams) 
{
    return atEnd(inStreams.stream);
}

//#ifdef __WITHCDR3THREADS__
//std::mutex MUTEX_readBlockOfRecords;
//#endif
//template<typename TSequencingSpec>
//long readBlockOfHighQualityRecords(
//        QueryDataCollection<TSequencingSpec> & queryDataCollection,
//        SeqInputStreams<TSequencingSpec> & inStreams,
//        CdrOptions const & options) 
//{
//    // ============================================================================
//    // Locking in case of multithreading
//    // ============================================================================
//#ifdef __WITHCDR3THREADS__
//    std::unique_lock<std::mutex> lock(MUTEX_readBlockOfRecords);
//#endif
//
//    unsigned nReads = 0;
//    StringSet<CharString> rejectIds;
//
//    if (inStreamsAtEnd(inStreams))
//        return -1;
//
//    while (!inStreamsAtEnd(inStreams))
//    {
//        FastqRecord<TSequencingSpec> fastqRecord;
//
//        try {
//            readRecord(fastqRecord, inStreams);
//        } catch (Exception const & e) {
//            std::cerr << "[ERR] Error reading FASTQ file: " << e.what() << std::endl;
//            return -1;
//        }
//
//        
//        // ============================================================================
//        // Truncate reads if requested so by the user
//        // ============================================================================
//
//        truncateRead(fastqRecord, options.trunkReads, options.vReadCrop);
//
//        // ============================================================================
//        // Ensure that the forward and reverse complement sequence have the same
//        // orientation.  Swap them if the -r option was specified. In case of
//        // single-end sequencing reverse complement the sequence if the -r option was
//        // specified.
//        // ============================================================================
//
//        syncOrientation(fastqRecord, options);
//
//        if (checkQuality(fastqRecord, options)) {
//            addRecord(queryDataCollection, fastqRecord);
//
//            ++nReads;
//        } else {
//            appendValue(rejectIds, fastqRecord.id);
//        }
//
//        if (nReads >= options.maxBlockSize)
//            break;
//
//    }
//    {
//#ifdef __WITHCDR3THREADS__
//        std::unique_lock<std::mutex> lock(MUTEX_ofstream_rejectlog);
//#endif
//        for (Iterator<StringSet<CharString> const, Rooted>::Type it = begin(rejectIds); !atEnd(it); goNext(it))
//            REJECTLOG << *it << '\t' << "AVERAGE_QUAL_FAIL" << std::endl;
//    }
//    return length(rejectIds);
//}



// Caching works only based on the VDJ read sequence right now! TODO 
TQueryDataSequence const getCacheStoreKey(QueryData<SingleEnd> const & qData, unsigned pos) 
{
    return qData.seqs[pos];
}

TQueryDataSequence const getCacheStoreKey(QueryData<PairedEnd> const & qData, unsigned pos) 
{
    QueryData<PairedEnd>::TSequence key;
    append(key, qData.revSeqs[pos]);
    append(key, qData.fwSeqs[pos]);
    return key;
}

#ifdef __WITHCDR3THREADS__
std::mutex MUTEX_processReads;
std::mutex PCR_ERR_STAT_MUTEX;
#endif
template <typename TQueryData, typename TGlobal>
String<AnalysisResult> analyseReads(
        TQueryData & queryData,       // [IN]  The query data with the reads to analyse
        TGlobal & global)             // [IN]  Global parameters and data
{
    // ============================================================================
    // Types
    // ============================================================================

    typedef typename TQueryData::TSequence                      TSequence;
    typedef String<AnalysisResult>                              TResults;
    typedef typename Infix<TSequence>::Type                     TInfix;
    typedef SegmentMatch<typename Infix<TSequence const>::Type> TSegmentMatch;

    // ============================================================================
    // Declarations
    // ============================================================================

    TQueryData newQueries;
    TResults results;
    String<bool> wasCached;
    resize(results, nRecords(queryData));
    resize(wasCached, nRecords(queryData), false);

    // ============================================================================
    // Find the best overlap alignments for the V and J segments
    // THREADS: Supposedly safe. Shares the reference sequences ({right,left}Segs)
    //          with other threads, however, should be read-only. Cannot be 
    //          declared const because the SWIFT implementation doesn't allow it.
    // ============================================================================

    StringSet<String<SegmentMatch<Infix<String<Dna5Q> const>::Type> > > leftMatches, rightMatches;

    // Find best matching V-segments

    findBestSegmentMatch(leftMatches, queryData, global, LeftOverlap());

    if (global.options.outputAligments)
    {
        for (unsigned x = 0; x < length(leftMatches); ++x)
        {
            std::cerr << "\n\n========= V ALIGNMENTS FOR READ " << x << " ===========\n\n";
            for (auto y : leftMatches[x])
            {
                int overlapLength = length(row(y.align,1));
                int nErrors = (overlapLength - y.score) / 2;
                double errRate = 1.0 * nErrors / overlapLength;
                std::cerr << getDescriptor(global.references.leftMeta[y.db]) << " (score=" << y.score << ", length=" << overlapLength << ", errRate=" << errRate << ")\n" << y.align;            
            }
        }
    }


    // Find best matching J-segments

    findBestSegmentMatch(rightMatches, queryData, global, RightOverlap());

    if (global.options.outputAligments)
    {
        for (unsigned x = 0; x < length(rightMatches); ++x)
        {
            std::cerr << "\n\n========= J ALIGNMENTS FOR READ " << x << " ===========\n\n";
            for (auto y : rightMatches[x])
            {
                std::cerr << getDescriptor(global.references.rightMeta[y.db]) << " (score=" << y.score << ")\n" << y.align;
            }
        }
    }

    // ============================================================================
    // For each newly calculated result the begin and end position of the cdr3
    // region is stored for caching purposes. To be able to address them by the
    // result (outer) id, the String has as many elements as there are
    // sequences. The same applies to the positions that were modified during
    // V/J reference overwrites.
    // ============================================================================

    String<std::pair<unsigned,unsigned> > cdrBeginEndPairs;
    resize(cdrBeginEndPairs, nRecords(queryData), std::make_pair(0u,0u));
    StringSet<String<unsigned> > modifiedPositions;
    resize(modifiedPositions, nRecords(queryData));

    // ============================================================================
    // Perform the analysis for each element that was not found in the cache
    // THREADS: Safe.
    // ============================================================================

    for (size_t i = 0; i<nRecords(queryData); ++i) {

        // ============================================================================
        // Reject the read if the V or J identification failed
        // ============================================================================

        if (length(leftMatches[i]) == 0 || length(rightMatches[i]) == 0) {
            results[i] = AnalysisResult(SEGMENT_MATCH_FAILED);
            continue;
        }

        TSequence modSeq = getVDJReadSequences(queryData)[i];

        // ============================================================================
        // Try to identify the CDR3 region and reject the read if that fails
        // ============================================================================

        TInfix cdrInfix(modSeq);
        {
            RejectReason reject = findCDR3Region(cdrInfix, leftMatches[i], rightMatches[i], global);
            if (reject) {
                results[i] = AnalysisResult(reject);
                continue;
            }
        }

        // ============================================================================
        // Reject the read if the boundaries are broken or if the Cys and the Phe are
        // in different reading frames
        // ============================================================================

        unsigned cdrBegin   = beginPosition(cdrInfix);
        unsigned cdrEnd     = endPosition(cdrInfix);

        cdrBeginEndPairs[i] = std::make_pair(cdrBegin, cdrEnd);

        if (cdrBegin >= cdrEnd) {
            results[i] = AnalysisResult(BROKEN_CDR_BOUNDARIES);
            continue;
        }

        if ((cdrEnd-cdrBegin)%3 != 0) {
            results[i] = AnalysisResult(OUT_OF_READING_FRAME);
            continue;
        }

        if ((cdrEnd-cdrBegin)/3 < global.options.minCDR3Length) {
            results[i] = AnalysisResult(CDR3_TOO_SHORT);
            continue;
        }

        // ============================================================================
        // Translate the CDR3 region and reject the read if any nonsense translation
        // appears in the aa sequence.
        // ============================================================================

        String<AminoAcid> aaCdr3;
        translate(aaCdr3, cdrInfix);
        bool nonsense = false;

        //TODO replace 'X' by '*'

        for (Iterator<String<AminoAcid>, Rooted>::Type it = begin(aaCdr3); !atEnd(it); goNext(it))
            if (*it == convert<AminoAcid>('X')) {
                results[i] = AnalysisResult(NONSENSE_IN_CDR3);
                nonsense = true;
                break;
            }

        if (nonsense)
            continue;


        // ============================================================================
        // We store the length and #errors of the alignments for the right /
        // left reference sequences - if it is identical for every optimally
        // matching reference sequence.
        // ============================================================================

        // Store the segment ids of the matching segments
        std::set<unsigned> vHits, jHits;
        for (typename Iterator<String<TSegmentMatch>, Rooted>::Type match = begin(leftMatches[i]); !atEnd(match); goNext(match)) {
            vHits.insert(match->db);
        }
        for (typename Iterator<String<TSegmentMatch>, Rooted>::Type match = begin(rightMatches[i]); !atEnd(match); goNext(match)) {
            jHits.insert(match->db);
        }

        // ============================================================================
        // If the V and J alignment error positions are unique, they are included in
        // the fully detailed output
        // ============================================================================

        String<int> leftErrPos, rightErrPos;
        unsigned leftOuterMatchLen = -1u;
        unsigned rightOuterMatchLen = -1u;
        bool leftUniqueErrPos = getErrorPositions(leftErrPos, leftOuterMatchLen, leftMatches[i], global.references.leftMeta, LeftOverlap());
        bool rightUniqueErrPos = getErrorPositions(rightErrPos, rightOuterMatchLen, rightMatches[i], global.references.rightMeta, RightOverlap());

        // ============================================================================
        // Assemble the tab separated string containing the detailed per read data
        // ============================================================================

        std::stringstream fullOut;
        if (global.outFiles._fullOutStream != NULL) {
            CharString leftString, rightString;
            _appendGeneList(leftString, vHits, global.references.leftMeta, global.options.mergeAllels);
            _appendGeneList(rightString, jHits, global.references.rightMeta, global.options.mergeAllels);
            String<AminoAcid> aaString;
            translate(aaString, cdrInfix);



            // CDR3 boundaries - the motif position is unique, as otherwise the read would have been rejected by now
            //                fullOut << "\t" << leftMatches[readId][0].motifPos << '\t' << rightMatches[readId][0].motifPos + 2;
            fullOut << "\t" << cdrBegin << '\t' << cdrEnd;

            // Add matching V segments
            fullOut << "\t" << leftString;

            // Add the V segment error positions 
            if (!leftUniqueErrPos) {
                fullOut << "\tNA\tNA";
            } else {
                fullOut << "\t";
                for (Iterator<String<int>, Rooted>::Type it = begin(leftErrPos); !atEnd(it); ) {
                    fullOut << *it;
                    goNext(it);
                    if (!atEnd(it))
                        fullOut << ",";
                    else
                        break;
                }
                fullOut << "\t" << leftOuterMatchLen;
            }


            // Add matching J segments
            fullOut << "\t" << rightString;

            // Add the J segment error positions 
            if (!rightUniqueErrPos) {
                fullOut << "\tNA\tNA";
            } else {
                fullOut << "\t";
                for (Iterator<String<int>, Rooted>::Type it = begin(rightErrPos); !atEnd(it); ) {
                    fullOut << *it;
                    goNext(it);
                    if (!atEnd(it))
                        fullOut << ",";
                    else
                        break;
                }
                fullOut << "\t" << rightOuterMatchLen;
            }

            fullOut << "\t" << cdrInfix << "\t" << aaString << std::endl;
        }

        // Store the result
        Clone<Dna5> c = { vHits, cdrInfix, jHits };
        String<uint64_t> cdrQualities;
        getQualityString(cdrQualities, cdrInfix);
        AnalysisResult cr(c, CharString(fullOut.str()), cdrQualities);

        results[i] = cr;

    }

    // ============================================================================
    // Write the detailed results for each read to a file if requested. Store the
    // new results in the cache store.
    // THREADS: *Not* safe. Alters the global clusterStore, writes to a global
    //          file.
    // ============================================================================

#ifdef __WITHCDR3THREADS__
    std::lock_guard<std::mutex> guard(MUTEX_processReads);
    std::lock_guard<std::mutex> guard2(MUTEX_ofstream_rejectlog);
    std::unique_lock<std::mutex> lock(CACHE_STORE_MUTEX);
#endif

    return results;
}

String<AnalysisResult> analyseReads(
        QueryDataCollection<SingleEnd> & qDatCol,       // [IN]  The query data collection with the reads to analyse
        CdrGlobalData<SingleEnd> & global)              // [IN]  Global parameters and data
{
    return analyseReads(qDatCol.queryData, global);
}

String<AnalysisResult> analyseReads(
        QueryDataCollection<PairedEnd> & qDatCol,       // [IN]  The query data collection with the reads to analyse
        CdrGlobalData<PairedEnd> & global)              // [IN]  Global parameters and data
{
    // Analyse the paired end data
    String<AnalysisResult> peRes = analyseReads(qDatCol.pairedQueryData, global);
    // Analyse the single end data
    String<AnalysisResult> seRes = analyseReads(qDatCol.singleQueryData, global);
    // Interlace correctly
    String<AnalysisResult> res;
    reserve(res, length(peRes)+length(seRes));
    size_t pe_idx = 0, se_idx = 0;
    Iterator<String<size_t> const, Rooted>::Type sePosIt = begin(qDatCol.sePositions);
    for (size_t i=0; i<nRecords(qDatCol); ++i)
    {
        if (sePosIt != end(qDatCol.sePositions) && *sePosIt == i)
        {
            ++sePosIt;
            appendValue(res, seRes[se_idx++]);
        } else {
            appendValue(res, peRes[pe_idx++]);
        }
    }
    return res;
}

#ifdef __WITHCDR3THREADS__
std::mutex MUTEX_take_from_multicollection;
std::mutex MUTEX_write_global_results;
#endif
template<typename TSeqSpec>
void processReads(
        String<AnalysisResult> & results,  // OUT: The map of clone counts to insert the results into
        ProgressBar& progBar,       //  IN: The ProgressBar object to report to
        FastqMultiRecordCollection<TSeqSpec> & collection,
        size_t & processed,
        CdrGlobalData<TSeqSpec> & global)           //  IN: The user specified parameters
{
    typedef QueryDataCollection<TSeqSpec>  TQueryDataCollection;
    while (true) { // Breaks when no more reads can be read from the input streams

        // ============================================================================
        // Extract block of records for analysis
        // ============================================================================

        String<FastqMultiRecord<TSeqSpec>*> todo;
        TQueryDataCollection qdataColl;
        size_t first_idx;
        { // Scope for lock
#ifdef __WITHCDR3THREADS__
            std::unique_lock<std::mutex> lock(MUTEX_take_from_multicollection);
#endif
            // Detect, how man items have to be processed
            size_t n_items = std::min<size_t>(length(collection.multiRecordPtrs) - processed, global.options.maxBlockSize);

            // Break loop if we are finished
            if (n_items == 0)
                break;

            // Store items that have to be processed
            first_idx = processed;
            for (size_t i = 0; i < n_items; ++i) {
                appendValue(todo, collection.multiRecordPtrs[i+processed]);
            }
            qdataColl = buildQDCollection(todo);
            processed += n_items;
        }

        // ============================================================================
        // Perform the actual analysis
        // ============================================================================

        String<AnalysisResult> results_block = analyseReads(qdataColl, global);
        progBar.updateAndPrint(length(todo));

        // ============================================================================
        // Process Analysis results
        // ============================================================================

        {
#ifdef __WITHCDR3THREADS__
            std::unique_lock<std::mutex> lock(MUTEX_write_global_results);
#endif
            for (AnalysisResult const & ar : results_block)
            {
                results[first_idx++] = ar;
            }
        }


    }
}

template <typename TSequencingSpec>
void splitAnalysisResults(
        std::map<Clone<Dna5>, ClusterResult> & nucCloneStore,
        String<RejectEvent> & rejectEvents,
        String<AnalysisResult> const & results,
        String<FastqMultiRecord<TSequencingSpec>*> const & recPtrs
        )
{
    nucCloneStore.clear();
    clear(rejectEvents);
    for (size_t i = 0; i<length(results); ++i)
    {
        AnalysisResult const & ar = results[i];
        FastqMultiRecord<TSequencingSpec> const & record = *(recPtrs[i]);
        if (!ar.reject) {
            // Increase the counter for this clone
            countNewClone(nucCloneStore, ar.clone, ar.cdrQualities);
        } else {
            for (CharString const & id : record.ids)
                appendValue(rejectEvents, RejectEvent(id, ar.reject));
        }
    }
}

template <typename TSequencingSpec>
void writeRDTFile(
        CdrGlobalData<TSequencingSpec> & global,
        String<AnalysisResult> const & results,
        String<FastqMultiRecord<TSequencingSpec>*> const & recPtrs
        )
{
    for (size_t i=0; i<length(results); ++i) {
        AnalysisResult const & ar = results[i];
        if (ar.reject)
            continue;
        FastqMultiRecord<TSequencingSpec> & rec = *recPtrs[i];
        for (CharString const & id : rec.ids) {
            *global.outFiles._fullOutStream << id << ar.fullOutSuffix;
        }
    }
}

template <typename TAlph>
size_t nClones(std::map<Clone<TAlph>, ClusterResult> const & cloneStore)
{
    typedef std::map<Clone<TAlph>, ClusterResult> TCloneStore;
    typedef typename TCloneStore::const_iterator  TIter;
    size_t cloneCount = 0;
    for (TIter ccount = cloneStore.begin(); ccount!=cloneStore.end(); ++ccount)
        cloneCount += ccount->second.count;
    return cloneCount;
}

template<typename TSequencingSpec>
String<AnalysisResult> runCDR3Analysis(
        FastqMultiRecordCollection<TSequencingSpec> & collection, //  IN: Input data
        CdrGlobalData<TSequencingSpec> & global                   //  IN: The user specified parameters
        )
{
    // ============================================================================
    // Status message
    // ============================================================================

    std::cerr << "===== Gene segment and CDR analysis" << std::endl;

    // ============================================================================
    // Initialize the status indicator within the command name if requested by the
    // user. Initialize the progress bar shown at the command prompt.
    // ============================================================================

    ProgressBar progBar(std::cerr, length(collection.multiRecordPtrs), 100, "      ");
    progBar.print_progress();

    // ============================================================================
    // Launch the analysis. processReads() reads a block of reads and analyses it.
    // If supported and requested by the user, multiple threads call
    // processReads(), the function and therefore each thread terminates when the
    // SequenceStream is exhausted, i.e. no more reads are available.
    // ============================================================================

    size_t processed = 0;
    String<AnalysisResult> results;
    resize(results, length(collection.multiRecordPtrs));
#ifdef __WITHCDR3THREADS__
    std::vector<std::thread> threads;
    for (int w=0; w < global.options.jobs; ++w)
        threads.push_back(std::thread(
                [&]() {
#endif
                processReads(results , progBar, collection, processed, global);
#ifdef __WITHCDR3THREADS__
                }
                ));
    for (std::thread & t : threads)
        t.join();
#endif
    progBar.clear();

    return results;
}

template <typename TStream>
void writeRejectLog(TStream & stream, String<RejectEvent> const & rejectEvents)
{
    for (RejectEvent const & rejectEvent : rejectEvents)
    {
        stream << rejectEvent.id << '\t' << _CDRREJECTS[rejectEvent.reason] << '\n';
    }
    stream << std::flush;
}

template<typename TSequencingSpec>
int runAnalysis(CdrGlobalData<TSequencingSpec> & global,
        String<RejectEvent> & rejectEvents,
        FastqMultiRecordCollection<TSequencingSpec> & collection)
{

    CdrOptions const & options = global.options;

    // ============================================================================
    // First pass over FASTQ input to determine read count and max length
    // ============================================================================

    std::ofstream *__aminoOutStream, *__nucOutStream;

    initOutFileStream(options.rlogPath, __rejectLog);
    initOutFileStream(options.aminoOut, __aminoOutStream);
    initOutFileStream(options.nucOut, __nucOutStream);
    initOutFileStream(options.fullOut, global.outFiles._fullOutStream);

    // ============================================================================
    // Write the CSV headers
    // ============================================================================

    PHREDSTAT << "length,errRate,end,ismatch";
    for (unsigned score = 0; score <= 60; ++score)
        PHREDSTAT << "," << score;
    PHREDSTAT << std::endl;

    POINTERSTREAM(global.outFiles._fullOutStream) << "seqId\tcdrBegin\tcdrEnd\tleftMatches\tleftErrPos\tleftMatchLen\trightMatches\trightErrPos\trightMatchLen\tcdrNucSeq\tcdrAASeq" << std::endl;

    // ============================================================================
    // Perform the analysis
    // ============================================================================

    std::clock_t beforeAnalysis = std::clock();

    String<AnalysisResult> results = runCDR3Analysis(collection, global);

    std::cerr << "  |-- Required cpu time: " << formatSeconds(double(std::clock() - beforeAnalysis) / CLOCKS_PER_SEC) << std::endl;;

    // ============================================================================
    // Fill Nucleotide Clone Store and reject events
    // ============================================================================

    typedef std::map<Clone<Dna5>, ClusterResult> TNucCloneStore;
    TNucCloneStore nucCloneStore;
    String<RejectEvent> newRejectEvents;
    splitAnalysisResults(nucCloneStore, newRejectEvents, results, collection.multiRecordPtrs);
    append(rejectEvents, newRejectEvents);

    // ============================================================================
    // Output info about how many of the reads were successfully analyzed
    // ============================================================================

    size_t cloneCount = nClones(nucCloneStore);
    std::string s = cloneCount == 1 ? "clone" : "clones";
    std::cerr << "  |-- " << cloneCount << " " << s << " could be identified." << std::endl;

    // ============================================================================
    // Write detailed per read output file if requested
    // ============================================================================

    if (global.outFiles._fullOutStream != NULL)
        writeRDTFile(global, results, collection.multiRecordPtrs);

    // ============================================================================
    // Write reject information if requested
    // ============================================================================

    if (__rejectLog != NULL)
        writeRejectLog(*__rejectLog, rejectEvents);

    // ============================================================================
    // Compute the average quality scores
    // ============================================================================

    std::cerr << "===== Computing average quality scores per position for each clonotype\n";
    for (TNucCloneStore::iterator it = nucCloneStore.begin(); it!=nucCloneStore.end(); ++it)
        it->second.buildAverageQScores();

    // ============================================================================
    // Drop allele information if requested so
    // ============================================================================

    if (global.options.mergeAllels) {
        std::cerr << "===== Removing allel information from clonotypes\n";
        unsigned oldClones = nucCloneStore.size();
        removeAllelInformation(nucCloneStore, global.references);
        std::cerr << "      Reduced from " << oldClones << " to " << nucCloneStore.size() << " clonotypes.\n";
    }

    // ============================================================================
    // Perform additional post processing steps if requested so
    // ============================================================================

    // =1= CLONOTYPE CLUSTERING

    runClonotypeClustering(global, nucCloneStore);

    // =2= POSTERIOR LOW QUALITY CLUSTER REMOVAL

    dropLowQualityClusters(nucCloneStore, options.qminclust);

    // =3= RESOLVE AMBIGUOUS V/J ASSIGNMENTS

    // Check if the user requested ambiguity resolution
    if (options.mergeIdenticalCDRs)
    {
        std::cerr << "===== Resolving V- and J-segment ambiguity (nucleotide-based)\n";
        resolveSegmentAmbiguity(nucCloneStore, global);
    }

    // ============================================================================
    // Write the output to the specified output files
    // ============================================================================

    if (__aminoOutStream != NULL) {

        std::map<Clone<AminoAcid>, ClusterResult> aaCloneStore;
        _translateClones(aaCloneStore, nucCloneStore);

        // Re-run the V/J ambiguity handling, since more clonotypes share the same
        // CDR3 sequence on the AA level
//        if (options.mergeIdenticalCDRs) {
//            mergeIdenticalCDRs(aaCloneStore);
//        }
        if (options.mergeIdenticalCDRs)
        {
            std::cerr << "===== Resolving V- and J-segment ambiguity (aminoacid-based)\n";
            resolveSegmentAmbiguity(aaCloneStore, global);
        }
        _writeClonotypeCounts(*__aminoOutStream, aaCloneStore, global);
    }

    if (__nucOutStream != NULL) {
        _writeClonotypeCounts(*__nucOutStream, nucCloneStore, global);
    }

    nucCloneStore.clear();
    closeOfStream(__aminoOutStream);
    closeOfStream(__nucOutStream);

    std::cerr << "===== Sorting output files\n";

    closeOfStream(__rejectLog);
    if (__rejectLog != NULL) {
        std::cerr << "  |-- Rejectlog\n";
        sortFile(options.rlogPath);
    }

    closeOfStream(global.outFiles._fullOutStream);
    if (global.outFiles._fullOutStream != NULL) {
        std::cerr << "  |-- Detailed output file\n";
        sortFile(options.fullOut, 1);
    }

    std::cerr << "===== All done. Terminating." << std::endl;

    return 0;
}

void autoTuneVSCFLength(CdrOptions & options, unsigned minReadLength) {
    minReadLength = options.trunkReads > 0 && options.trunkReads < minReadLength ? options.trunkReads : minReadLength;
    if (options.vSCFLengthAuto) {
        if (minReadLength >= 120) options.vSCFLength = 20;
        else if (minReadLength < 100) options.vSCFLength = 10;
        else options.vSCFLength = 10 + (10-((120 - minReadLength)/2));
    }
}


template <typename TStream>
void printStats(TStream & stream, String<RejectEvent> const & rejectEvents,
        InputInformation const & ii,
        BarcodeStats const & stats)
{
    stream <<
        "  |   ...... Number of reads read: " << ii.totalReadCount << '\n' <<
        "  |   .................. Rejected: " << length(rejectEvents) << '\n' <<
        "  |   .. Max read length (passed): " << ii.maxReadLength << '\n' << 
        "  |   .. Min read length (passed): " << ii.minReadLength << '\n' <<
        "  |   Number of unique read pairs: " << stats.nTotalUniqueReads << '\n';
}

inline void printPerBarcodeStats(BarcodeStats const & stats)
{
    std::cout << "Barcode\tnReads\tnUniqueReads\n";
    for (unsigned i=0; i<length(stats.bcSeqs); ++i)
        std::cout << stats.bcSeqs[i] << "\t" << stats.nReads[i] << "\t" << stats.nUniqueReads[i] << std::endl;
    std::cerr << "      " << stats.nTotalReads << " read pairs (" << stats.nTotalUniqueReads << " unique)" << std::endl;
}

template <typename TSequencingSpec>
int main_generic(CdrGlobalData<TSequencingSpec> & global, CdrOptions & options, CdrReferences & references) {

    ignoreUnusedVariableWarning(references);

    // ============================================================================
    // READ FASTQ FILES
    // ============================================================================

    std::cerr << "===== PROCESSING INPUT FILES\n";
    // Target data structures
    InputInformation inputInformation;
    FastqMultiRecordCollection<TSequencingSpec> collection;
    String<RejectEvent> rejectEvents;
    // Read data
    readRecords(collection, inputInformation, rejectEvents, global.input, global.options);
    // Print information
    BarcodeStats stats = getBarcodeStats(collection);
    printStats(std::cerr, rejectEvents, inputInformation, stats);

    // ============================================================================
    // TUNE V-SCF LENGTH, READ REFERENCE SEGMENT SEQUENCES
    // ============================================================================

    autoTuneVSCFLength(options, inputInformation.minReadLength);
    readAndPreprocessReferences(references, global.options);

    std::cerr << "  |-- Read " << length(global.references.leftSegs) << " reference V segments." << std::endl;
    std::cerr << "  |-- Read " << length(global.references.rightSegs) << " reference J segments." << std::endl;


    // ============================================================================
    // COLLAPSE SIMILAR BARCODES
    // ============================================================================

    std::cerr << "===== BARCODE CORRECTION\n";
    std::cerr << "  |-- Joining similar barcodes (" << options.barcodeMaxError << " errors allowed)" << std::endl;

    clusterBarcodeSequences(collection, options);
    stats = getBarcodeStats(collection);
    std::cerr <<
        "  |   ...... Number of reads read: " << stats.nTotalReads << '\n' <<
        "  |   Number of unique read pairs: " << stats.nTotalUniqueReads << '\n';
    compact(collection);

    // ============================================================================
    // CORRECT READS BASED ON BARCODES
    // ============================================================================

    std::cerr << "  |-- Performing barcode correction" << std::endl;
    barcodeCorrection(collection, options);
    stats = getBarcodeStats(collection);
    std::cerr <<
        "  |   ...... Number of reads read: " << stats.nTotalReads << '\n' <<
        "  |   Number of unique read pairs: " << stats.nTotalUniqueReads << '\n';
    compact(collection);

    // ============================================================================
    // DROP BARCODE INFORMATION
    // ============================================================================

    FastqMultiRecordCollection<TSequencingSpec> noBcCollection;

    for (FastqMultiRecord<TSequencingSpec> * const recPtr : collection.multiRecordPtrs)
    {
        if (recPtr == NULL || recPtr->ids.empty())
            continue;
        FastqMultiRecord<TSequencingSpec> recCopy = *recPtr;
        recCopy.bcSeq = "";
        mergeRecord(noBcCollection, recCopy);
    }

    std::cerr << "  |-- Dropping barcode information" << std::endl;
    stats = getBarcodeStats(noBcCollection);
    std::cerr <<
        "      ...... Number of reads read: " << stats.nTotalReads << '\n' <<
        "      Number of unique read pairs: " << stats.nTotalUniqueReads << '\n';

    // ============================================================================
    // RUN THE CLONOTYING ANALYSIS
    // ============================================================================

    return runAnalysis(global, rejectEvents, noBcCollection);
}

#endif  // #ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_CDR3FINDER_H_

/* vim: set sw=4 sts=4 ts=8 spell spelllang=en expandtab: */
