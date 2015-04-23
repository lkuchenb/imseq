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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_COLLECTION_UTILS_H_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_COLLECTION_UTILS_H_

#include <map>

template<typename TC1, typename TC2>
bool shareElement(TC1 const & c1, TC2 const & c2) {
    for (typename TC1::const_iterator it1 = c1.begin(); it1!=c1.end(); ++it1)
        for (typename TC2::const_iterator it2 = c2.begin(); it2!=c2.end(); ++it2)
            if (*it1 == *it2)
                return true;
    return false;
}

template<typename TKey, typename TValue>
bool shareKey(std::map<TKey,TValue> const & m1, std::map<TKey,TValue> const & m2) {
    for (typename std::map<TKey,TValue>::const_iterator it1 = m1.begin(); it1!=m1.end(); ++it1)
        for (typename std::map<TKey,TValue>::const_iterator it2 = m2.begin(); it2!=m2.end(); ++it2)
            if (it1->first == it2->first)
                return true;
    return false;
}

#endif  // #ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_PROGRESS_BAR_H_
