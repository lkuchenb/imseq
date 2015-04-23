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
        {'K','N','K','N','X'},
        {'T','T','T','T','X'},
        {'R','S','R','S','X'},
        {'I','I','M','I','X'},
        {'X','X','X','X','X'}
    },
    {
        {'Q','H','Q','H','X'},
        {'P','P','P','P','X'},
        {'R','R','R','R','X'},
        {'L','L','L','L','X'},
        {'X','X','X','X','X'}
    },
    {

        {'E','D','E','D','X'},
        {'A','A','A','A','X'},
        {'G','G','G','G','X'},
        {'V','V','V','V','X'},
        {'X','X','X','X','X'}
    },
    {

        {'X','Y','X','Y','X'},
        {'S','S','S','S','X'},
        {'X','C','W','C','X'},
        {'L','F','L','F','X'},
        {'X','X','X','X','X'}
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
 * Nonsense-codons as well as triplets containing 'N's are translated with 'X'.
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

#endif  // #ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_AA_TRANSLATE_
