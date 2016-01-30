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

#ifndef CDR3FINDER_FILE_UTILS
#define CDR3FINDER_FILE_UTILS

#include <seqan/sequence.h>
#include <vector>
#include <iostream>
#include <algorithm>

inline void sortFile(CharString path, unsigned keepLines = 0) {
    std::vector<std::string> fileContent, fileContentKeep;
    std::ifstream ifs(toCString(path));
    std::string s;
    while (keepLines > 0 && std::getline(ifs, s)) {
        --keepLines;
        fileContentKeep.push_back(s);
    }
    while (std::getline(ifs, s)) {
        fileContent.push_back(s);
    }
    std::sort(fileContent.begin(), fileContent.end());
    ifs.close();
    std::ofstream ofs(toCString(path));
    for (std::vector<std::string>::const_iterator it = fileContentKeep.begin(); it!=fileContentKeep.end(); ++it)
        ofs << *it << '\n';
    for (std::vector<std::string>::const_iterator it = fileContent.begin(); it!=fileContent.end(); ++it)
        ofs << *it << '\n';
    ofs.close();
}

#endif
