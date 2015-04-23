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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_CDR_UTILS_H_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_CDR_UTILS_H_

#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <seqan/sequence.h>

template<typename T>
T seqanStringSum(seqan::String<T> const & str) {
    T sum = 0;
    for (typename seqan::Iterator<seqan::String<T> const,seqan::Rooted>::Type it = begin(str); !atEnd(it); goNext(it))
        sum += *it;
    return sum;
}

template<typename TKey, typename TValue>
std::vector<TKey> mapKeys(std::map<TKey, TValue> const & map) {
    typedef typename std::map<TKey, TValue>     TMap;
    std::vector<TKey> res;
    for (typename TMap::const_iterator it = map.begin(); it!=map.end(); ++it)
        res.push_back(it->first);
    return res;
}

template<typename TKey, typename TValue>
std::vector<TKey> mapKeyVector(std::map<TKey, TValue> const & map) 
{
    typedef typename std::map<TKey, TValue>     TMap;
    std::vector<TKey> vec;
    for (typename TMap::const_iterator it = map.begin(); it!=map.end(); ++it)
        vec.push_back(it->first);
    return vec;
}


template<typename T>
bool collectionOverlapSTL(T const & a, T const & b) {
    for (typename T::const_iterator aIt = a.begin(); aIt!=a.end(); ++aIt)
        for (typename T::const_iterator bIt = b.begin(); bIt!=b.end(); ++bIt)
            if (*aIt == *bIt) return true;
    return false;
}

template<typename TSTLContainer>
std::string toSeparatedStringSTL(TSTLContainer const & original, char beginChar='{', char endChar = '}', char sepChar=',') {

    std::stringstream res;

    res << beginChar;
    if (original.size() > 0)
        for (typename TSTLContainer::const_iterator it = original.begin() ; ; ) {
            res << *it;
            ++it;
            if (it==original.end())
                break;
            else
                res << sepChar;
        }
    res << endChar;


    return(res.str());
}

template<typename T>
std::set<T> setIntersect(std::set<T> const & a, std::set<T> const & b) {
    std::set<T> res;
    for (typename std::set<T>::const_iterator aIt = a.begin(); aIt!=a.end(); ++aIt)
        if (b.find(*aIt) != b.end())
            res.insert(*aIt);
    return res;
}

template<typename T>
bool doIntersect(T const & a, T const & b) {
    for (typename seqan::Iterator<T const>::Type it = seqan::begin(a); it!=seqan::end(a); ++it)
        for (typename seqan::Iterator<T const>::Type iit = seqan::begin(b); iit!=seqan::end(b); ++iit)
            if (*iit == *it)
                return true;
    return false;
}

template<typename T>
bool doIntersectSTL(T const & a, T const & b) {
    for (typename T::const_iterator it = a.begin(); it!=a.end(); ++it)
        for (typename T::const_iterator iit = b.begin(); iit!=b.end(); ++iit)
            if (*iit == *it)
                return true;
    return false;
}

template<typename TKey, typename TValue>
TValue& mapMustHaveKey(std::map<TKey, TValue> & m, TKey const & key, char const * msg = NULL, int returnCode = 1) {
    typename std::map<TKey, TValue>::iterator it = m.find(key);
    if (it == m.end()) {
        std::string msgStr = msg == NULL ? "" : std::string(msg);
        std::cerr << '\n' << msgStr << '\n';
        exit(returnCode);
    }
    return it->second;
}

template<typename TKey, typename TValue>
TValue const & mapMustHaveKey(std::map<TKey, TValue> const & m, TKey const & key, char const * msg = NULL, int returnCode = 1) {
    typename std::map<TKey, TValue>::const_iterator it = m.find(key);
    if (it == m.end()) {
        std::string msgStr = msg == NULL ? "" : std::string(msg);
        std::cerr << '\n' << msgStr << " ### PLEASE REPORT THIS ERROR\n";
        exit(returnCode);
    }
    return it->second;
}

inline void exitWithError(char const * msg = NULL, int returnCode = 1) {
    if (msg != NULL)
        std::cerr << '\n' << msg << " ### PLEASE REPORT THIS ERROR\n";
    exit(returnCode);
}
#endif
