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

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_AA_TRANSLATE_
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_AA_TRANSLATE_

#include <cctype>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*
 * 3-dim array representing the genetic code, ordered by the integer
 * representation of the Dna alphabet
 */
static const AminoAcid _geneticCode[5][5][5] = {
    {
        {'K','N','K','N','*'},
        {'T','T','T','T','*'},
        {'R','S','R','S','*'},
        {'I','I','M','I','*'},
        {'*','*','*','*','*'}
    },
    {
        {'Q','H','Q','H','*'},
        {'P','P','P','P','*'},
        {'R','R','R','R','*'},
        {'L','L','L','L','*'},
        {'*','*','*','*','*'}
    },
    {

        {'E','D','E','D','*'},
        {'A','A','A','A','*'},
        {'G','G','G','G','*'},
        {'V','V','V','V','*'},
        {'*','*','*','*','*'}
    },
    {

        {'*','Y','*','Y','*'},
        {'S','S','S','S','*'},
        {'*','C','W','C','*'},
        {'L','F','L','F','*'},
        {'*','*','*','*','*'}
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/*
 * Translates a DNA or RNA String to an amino acid string. The translation
 * is performed on the first reading frame, trailing nucleotides are ignored.
 * Nonsense-codons as well as triplets containing 'N's are translated with '*'.
 */
template<typename TSpecIn, typename TSource>
void translate(String<AminoAcid, TSpecIn>& target, const TSource& nucSeq, unsigned short rFrame = 0) {
    clear(target);
    rFrame = rFrame % 3;
    if (length(nucSeq) < 3)
        return;
    reserve(target,length(nucSeq)/3);
    unsigned nucLength = length(nucSeq);
    nucLength -= (nucLength % 3);
    for (unsigned triplet = 0; triplet< nucLength-2; triplet+=3)
        appendValue(target, _geneticCode[(unsigned short)nucSeq[triplet]][(unsigned short)nucSeq[triplet+1]][(unsigned short)nucSeq[triplet+2]]);
}

/*
 * Translates a DNA or RNA String to an amino acid string. The translation
 * is performed on the first reading frame, trailing nucleotides are ignored.
 * Nonsense-codons as well as triplets containing 'N's are translated with '*'.
 * This version takes into account the case of the input string, i.e. it works
 * on char strings and outputs a lower case amino acid if the codon contained
 * at least one lower case base.
 */
template<typename TSpecOut, typename TSpecIn>
void translateCase(String<char, TSpecOut> & target, const String<char, TSpecIn> & nucSeq, unsigned short rFrame = 0) {
    clear(target);
    rFrame = rFrame % 3;
    if (length(nucSeq) < 3)
        return;
    reserve(target,length(nucSeq)/3);
    unsigned nucLength = length(nucSeq);
    nucLength -= (nucLength % 3);
    for (unsigned triplet = 0; triplet< nucLength-2; triplet+=3)
    {
        AminoAcid aa = _geneticCode[(unsigned short)convert<Dna5>(nucSeq[triplet])][(unsigned short)convert<Dna5>(nucSeq[triplet+1])][(unsigned short)convert<Dna5>(nucSeq[triplet+2])];
        bool lower = false;
        for (unsigned short i = 0; i < 3; ++i)
            if (std::islower(nucSeq[triplet+i]))
            {
                lower = true;
                break;
            }
        char c = convert<char>(aa);
        appendValue(target, lower ? std::tolower(c) : std::toupper(c));
    }
}

#endif  // #ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_AA_TRANSLATE_
