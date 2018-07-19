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

// ============================================================================
// HEADER FILE DESCRIPTION
// ============================================================================
// In this header file the functions for matching V and J segment sequences
// against the input sequences (reads) are defined.
// ============================================================================

#ifndef SANDBOX_LKUCHENB_APPS_CDR3FINDER_VJ_MATCHING_H
#define SANDBOX_LKUCHENB_APPS_CDR3FINDER_VJ_MATCHING_H

#include <cstdlib>
#include <climits>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/align_extend.h>
#include "cdr3_cli.h"
#include "globalData.h"
#include "overlap_specs.h"
#include "extdir_oldir_conversion.h"

/********************************************************************************
 * STRUCTS AND CLASSES
 *******************************************************************************/

template<typename TSeq>
struct RefMatch {
    typedef Align<TSeq> TAlign;
    int                 score;
    TAlign              align;

    RefMatch(int score_, Align<TSeq> align_) : 
        score(score_), align(align_) { }

    RefMatch() : score(0) { }

};

template<typename TSeq>
struct SegmentMatch : RefMatch<TSeq> {
    unsigned db;

    SegmentMatch(unsigned db_, int score_, Align<TSeq> align_) : RefMatch<TSeq>(score_, align_), db(db_) { }
};

struct CandidateCoreSegmentMatch
{
    unsigned coreSegId;
    unsigned readBeginPos;
    unsigned readEndPos;

    CandidateCoreSegmentMatch(unsigned coreSegId, unsigned readBeginPos, unsigned readEndPos)
        : coreSegId(coreSegId), readBeginPos(readBeginPos), readEndPos(readEndPos) {}
};

/********************************************************************************
 * FUNCTIONS
 *******************************************************************************/

/**
 * Constructs an ordered container with dependant infix segments.
 * The infixes have length 'length' and start at position 'start' relative to
 * the motif position.
 */
const unsigned BUILD_SEGMENT_CORE_FRAGMENTS_GOOD = -1u;
template<typename TStringSet>
unsigned buildSegmentCoreFragments(
        TStringSet & scfSet,                                    // [OUT] Set of segment core fragments
        StringSet<String<unsigned> > & scfToSegId,              // [OUT] Set of segment ids for each core fragment
        String<unsigned> & segToScfId,                          // [OUT] Mapping of segment IDs to SCF ids
        String<BeginEndPos<unsigned> > & scfBeginEndPos,        // [OUT] The SCF begin and end positions for each segment
        TStringSet const & segSeqs,                             // [IN]  The segment sequences
        String<SegmentMeta> const & segMeta,                    // [IN]  Meta information for the segments
        int const start,                                        // [IN]  The start position relative to the motif position
        unsigned const infixLength)                             // [IN]  The length of the core fragment
{
    typedef typename Value<TStringSet>::Type            TSequence;
    typedef typename Infix<TSequence const>::Type        TInfix;
    typedef std::map<TSequence, String<unsigned> >        TResMap;

    SEQAN_CHECK(length(segSeqs)==length(segMeta), "ERROR 1012 - length mismatch in buildSegmentCoreFragments(). Please report this error.");

    TResMap resMap;
    clear(scfBeginEndPos);
    for (unsigned i=0; i<length(segSeqs); ++i) {
        int beginPos = segMeta[i].motifPos+start;
        int endPos = beginPos + infixLength;

        if (beginPos >= 0 && beginPos < static_cast<long>(length(segSeqs[i])) && endPos > beginPos && endPos <= static_cast<long>(length(segSeqs[i]))) {
            TSequence sequenceInfix = TInfix(segSeqs[i], beginPos, endPos);
            appendValue(resMap[sequenceInfix], i);
            appendValue(scfBeginEndPos, BeginEndPos<unsigned>(beginPos, endPos));
        } else {
            clear(scfSet);
            return i;
        }
    }
    clear(scfSet);
    clear(scfToSegId);
    clear(segToScfId);
    resize(segToScfId, length(segSeqs), -1u);
    reserve(scfSet, resMap.size());
    reserve(scfToSegId, resMap.size());
    for (typename TResMap::const_iterator resMapIt = resMap.begin(); resMapIt!=resMap.end(); ++resMapIt) {
        appendValue(scfSet, resMapIt->first);
        appendValue(scfToSegId, resMapIt->second);
        for (unsigned segId : resMapIt->second)
            segToScfId[segId] = length(scfSet) - 1;
    }
    return BUILD_SEGMENT_CORE_FRAGMENTS_GOOD;
}

template<typename TAlign>
void clipOuterGaps(TAlign & align) {
    typedef typename Row<TAlign>::Type          TRow;

    TRow & readRow = row(align, 0);
    TRow & segRow = row(align, 1);

    // Clip non-matching head and tail
    int firstMatch = 0;
    while (firstMatch < static_cast<int>(length(readRow)) && (isGap(readRow, firstMatch) || isGap(segRow, firstMatch) || segRow[firstMatch]!=readRow[firstMatch]))
        ++firstMatch;

    int lastMatch = length(readRow)-1;
    while (lastMatch >= 0 && (isGap(readRow, lastMatch) || isGap(segRow, lastMatch) || segRow[lastMatch]!=readRow[lastMatch]))
        --lastMatch;

    // If there is at least one match, do the clipping
    if (firstMatch <= lastMatch)
    {
        setClippedEndPosition(readRow, clippedEndPosition(readRow) - (static_cast<int>(length(readRow)) - lastMatch - 1));
        setClippedEndPosition(segRow, clippedEndPosition(segRow) - (static_cast<int>(length(segRow)) - lastMatch - 1));
        setClippedBeginPosition(readRow, clippedBeginPosition(readRow) + firstMatch);
        setClippedBeginPosition(segRow, clippedBeginPosition(segRow) + firstMatch);
    }
}

/**
 * Clips a pairwise alignment such that there are no leading or trailing gaps
 * and computes the positions tuple required for extendAlignment().
 */
template<typename TInteger, typename TAlign, typename TOverlapDirection>
void clipAndComputePositions(Tuple<TInteger, 4> & positions, TAlign & align, TOverlapDirection const) {
    typedef typename Row<TAlign>::Type          TRow;

    TRow & readRow = row(align, 0);
    TRow & segRow = row(align, 1);

    // Clip non-matching head and tail
    int firstMatch = 0;
    while (firstMatch < static_cast<int>(length(readRow)) && (isGap(readRow, firstMatch) || isGap(segRow, firstMatch) || segRow[firstMatch]!=readRow[firstMatch]))
        ++firstMatch;

    int lastMatch = length(readRow)-1;
    while (lastMatch >= 0 && (isGap(readRow, lastMatch) || isGap(segRow, lastMatch) || segRow[lastMatch]!=readRow[lastMatch]))
        --lastMatch;

    // If there is at least one match, do the clipping
    if (firstMatch <= lastMatch)
    {
        setClippedEndPosition(readRow, clippedEndPosition(readRow) - (static_cast<int>(length(readRow)) - lastMatch - 1));
        setClippedEndPosition(segRow, clippedEndPosition(segRow) - (static_cast<int>(length(segRow)) - lastMatch - 1));
        setClippedBeginPosition(readRow, clippedBeginPosition(readRow) + firstMatch);
        setClippedBeginPosition(segRow, clippedBeginPosition(segRow) + firstMatch);
    }
        
    // Compute the begin and end positions in the hosts
    TInteger readBegin  = toSourcePosition(readRow, 0) + beginPosition(source(readRow));
    TInteger readEnd    = endPosition(readRow) + beginPosition(source(readRow));
    TInteger segBegin   = toSourcePosition(segRow, 0) + beginPosition(source(segRow));
    TInteger segEnd     = endPosition(segRow) + beginPosition(source(segRow));

    // Store the positions
    positions[0] = readBegin;
    positions[1] = segBegin;
    positions[2] = readEnd;
    positions[3] = segEnd;
}


/*
 * Computes the best score for a semi-global alignment with Levenshtein distance, given a text, a
 * pattern and a maximum number of errors. The score is returned and the end position of the best
 * match is written to 'endPos'. BANDED VERSION
 */
template<typename TText, typename TPattern>
int semiglobalAlignmentScoreMyersBanded(int & endPos, TText & text, TPattern & pattern, int maxErrors) {
    typedef Myers<AlignTextBanded<FindInfix, NMatchesNone_, NMatchesNone_>, True, void> TAlgorithmSpec;
    typedef PatternState_<TPattern, TAlgorithmSpec> TPatternState;
    typedef Finder<TText> TFinder;

    TFinder finder(text);
    TPatternState patternState;
    int minErrors = MaxValue<int>::VALUE;

    endPos = 0;

    while (find(finder, pattern, patternState, -maxErrors))
    {
        typename Position<TFinder>::Type currentEnd = position(finder) + 1;
        int currentErrors = -getScore(patternState);

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = currentEnd;

        }
    }
    return -minErrors;
}


/*
 * Computes the best score for a semi-global alignment with Levenshtein distance, given a text, a
 * pattern and a maximum number of errors. The score is returned and the end position of the best
 * match is written to 'endPos'. UNBANDED VERSION
 */
template<typename TText, typename TPattern>
int semiglobalAlignmentScoreMyers(int & endPos, TText & text, TPattern & pattern, int maxErrors) {
    typedef Pattern<CharString, Myers<> > TPatternState;
    typedef Finder<TText> TFinder;

    TFinder finder(text);
    TPatternState patternState(pattern);
    int minErrors = MaxValue<int>::VALUE;

    endPos = 0;

    while (find(finder, patternState, -maxErrors))
    {
        typename Position<TFinder>::Type currentEnd = position(finder) + 1;
        int currentErrors = -getScore(patternState);

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = currentEnd;

        }
    }
    return -minErrors;
}

template <typename TSequence>
String<CandidateCoreSegmentMatch> verifySCFHits(
        TSequence const & readSeq,
        StringSet<TSequence> const & scfSequences,
        std::map<unsigned, BeginEndPos<long long> > scfInfixPositions,
        int const maxErrors)
{
    typedef std::map<unsigned, BeginEndPos<long long> > TInfixPosMap;
    typedef typename Infix<TSequence const>::Type       TInfix;

    String<CandidateCoreSegmentMatch> ccsms;

    typedef TInfixPosMap::mapped_type TBeginEndPos;

    // Find all positions where the SCF matches within the constraints
    // specified by the user (# of errors)
    for (TInfixPosMap::const_iterator scfInfixPosition = scfInfixPositions.begin(); scfInfixPosition != scfInfixPositions.end(); ++scfInfixPosition) {
        unsigned coreSegId = scfInfixPosition->first;
        TBeginEndPos const & infixPos = scfInfixPosition->second;
        // Build the infix-SCF alignment
        TInfix readInfix(readSeq, infixPos.beginPos, infixPos.endPos);
        // "Pseudoinfix" to maintain type equality
        TInfix segInfix(scfSequences[coreSegId]);

        // Use Myers finder to find all locations within constraints
        Pattern<TInfix, Myers<> > myersPattern(segInfix);
        Finder<TInfix> finder(readInfix);

        int bestScore = MinValue<int>::VALUE;
        String<TBeginEndPos> toExpandInfixes;
        while (find(finder, myersPattern, -maxErrors))
        {
            if (getScore(myersPattern) < bestScore)
                continue;
            if (getScore(myersPattern) > bestScore)
            {
                bestScore = getScore(myersPattern);
                clear(toExpandInfixes);
            }
            TBeginEndPos bePos(MaxValue<long long>::VALUE, endPosition(finder));
            while (findBegin(finder, myersPattern, getScore(myersPattern)))
                if (bePos.beginPos > static_cast<long long>(beginPosition(finder)))
                    bePos.beginPos = beginPosition(finder);
            appendValue(toExpandInfixes, bePos);
        }

        // Store candidate match
        for (Iterator<String<TBeginEndPos>,Rooted>::Type it = begin(toExpandInfixes); !atEnd(it); goNext(it))
        {
            CandidateCoreSegmentMatch ccsm(coreSegId, it->beginPos + infixPos.beginPos, it->endPos + infixPos.beginPos);
            appendValue(ccsms, ccsm);
        }
    }
    return ccsms;
}

template <typename TSequence, typename TSegmentPattern>
std::map<unsigned, BeginEndPos<long long> > filterSCFs(
        TSequence & readSeq,
        TSegmentPattern & segmentPattern,
        double const errRate,
        std::set<unsigned> const * scfIdLim = nullptr)
{
    typedef std::map<unsigned, BeginEndPos<long long> >                 TInfixPosMap;
    typedef Finder<TSequence, Swift<SwiftSemiGlobal> >                  TFinder;

    TInfixPosMap scfInfixPositions;
    TFinder readFinder(readSeq);

    // No legal SCFS? Return empty result
    if (scfIdLim != nullptr && scfIdLim->empty())
        return scfInfixPositions;

    // Identify the windows within the read we have to investigate for potential
    // core fragment matches
    while (find(readFinder, segmentPattern, errRate, 0)) {
        unsigned coreSegId  = readFinder.curHit->ndlSeqNo;

        // Skip if a limited set was specified and the id is not in it
        if (scfIdLim != nullptr && scfIdLim->find(coreSegId) == scfIdLim->end())
            continue;

        long long beginPos  = readFinder.curHit->hstkPos;
        long long endPos    = beginPos + readFinder.curHit->bucketWidth;
        // 0 is the lower bound for our diagonal, since we want a semi global alignment
        beginPos            = beginPos < 0 ? 0 : beginPos;
        endPos              = endPos < 0 ? 0 : endPos;
        // Length of readSeq is the upper bound of our diagonal, since we want a semi global alignment
        beginPos            = beginPos < (long long)length(readSeq) ? beginPos : (long long)length(readSeq) - 1;
        endPos              = endPos <= (long long)length(readSeq) ? endPos : (long long)length(readSeq);
        // Check if there was a previous hit on the same segment
        TInfixPosMap::iterator scfInfixPosition = scfInfixPositions.find(coreSegId);
        if (scfInfixPosition != scfInfixPositions.end()) {
            BeginEndPos<long long> & infixPos       = scfInfixPosition->second;
            infixPos.beginPos = infixPos.beginPos > beginPos ? beginPos : infixPos.beginPos;
            infixPos.endPos   = infixPos.endPos < endPos ? endPos : infixPos.endPos;
        } else {
            BeginEndPos<long long> infixPos(beginPos, endPos);
            if (infixPos.endPos > infixPos.beginPos)
                scfInfixPositions[coreSegId] = infixPos;
        }
    }

    return scfInfixPositions;
}


template <typename TStringSet>
void findCandidateCoreSegments(
        StringSet<String<CandidateCoreSegmentMatch> > & results,        // [OUT] One String of candidate SCFs for each input read sequence
        TStringSet const & readSequences,                               // [IN]  The input read sequences
        TStringSet const & scfSequences,                                // [IN]  The core fragments of the segment sequences
        int const maxErrors                                             // [IN]  Maximum #errors for the core fragment alignment
        )
{
    // Extract types from passed arguments
    typedef typename Value<TStringSet>::Type                            TSequence;
    typedef typename Value<TSequence>::Type                             TAlphabet;

    // Types required for filtering
    typedef Shape<TAlphabet, SimpleShape>                               TShape;
    typedef Index<TStringSet, IndexQGram<TShape, OpenAddressing> >      TIndex;
    typedef Pattern<TIndex, Swift<SwiftSemiGlobal> >                    TPattern;
    typedef typename Infix<TSequence const>::Type                       TInfix;
    typedef Align<TInfix>                                               TAlign;
    typedef std::map<unsigned, BeginEndPos<long long> >                 TInfixPosMap;

    // Constants
    Score<int,Simple> const sc(1,-1,-1);

    // Pattern is constructed over segment core fragments, Finder is called on read sequences
    double const errRate = 1.0 * maxErrors / length(scfSequences[0]); // The length of all SCFs is the same
    TIndex segmentIndex(scfSequences);
    resize(indexShape(segmentIndex), length(scfSequences[0]) / (maxErrors+1));
    TPattern segmentPattern(segmentIndex);

    clear(results);
    reserve(results, length(readSequences));

    // Iterate through the read sequences and find potential core fragment matches
    for (typename Iterator<TStringSet const, Rooted>::Type readIt = begin(readSequences); !atEnd(readIt); goNext(readIt)) {
        // Need to give up const here. Maybe this won't be necessary some time in the 
        // distant future... 
        TSequence & readSeq = const_cast<TSequence &>(*readIt);

        // Filter
        TInfixPosMap scfInfixPositions = filterSCFs(readSeq, segmentPattern, errRate);

        // Verification
        String<CandidateCoreSegmentMatch> ccsms = verifySCFHits(readSeq, scfSequences, scfInfixPositions, maxErrors);

        appendValue(results, ccsms);
    }
}

template <typename TRow>
void _refineTerminalGaps(TRow & rowA, TRow & rowB)
{
    typedef typename Position<TRow>::Type TRowPos;

    // Left
    TRowPos vPos    = 1 + toViewPosition(rowB, length(source(rowB)) - 1);
    TRowPos maxGaps = length(rowB) - vPos;
    unsigned gaps = 0;
    for (; vPos > 0 && isGap(rowA, vPos - 1); --vPos)
        ++gaps;
    gaps = gaps > maxGaps ? maxGaps : gaps;
    if (gaps>0)
        removeGaps(rowA, vPos, gaps);

    // Right
    vPos = toViewPosition(rowA, 0);
    maxGaps = vPos;
    gaps = 0;
    for (; vPos < length(rowB) && isGap(rowB, vPos); ++vPos)
        ++gaps;
    gaps = gaps > maxGaps ? maxGaps : gaps;
    if (gaps>0)
    {
        removeGaps(rowB, toViewPosition(rowA, 0), gaps);
        removeGaps(rowA, 0, gaps);
    }
}

/**
 * Remove score indifferent gaps at the end of segment - read alignments
 *
 * This function takes care of the following type of scenario:
 *
 * *************-********       *********************
 * |||||||||||||           ==>  |||||||||||||
 * **************--------       **************-------
 *
 * which is score indifferent due to free end gaps.
 */
template <typename TAlign>
void refineTerminalGaps(TAlign & align)
{
    _refineTerminalGaps(row(align, 0), row(align, 1));
    _refineTerminalGaps(row(align, 1), row(align, 0));
}

template <typename TAlign, typename TReadSequence>
int extendToOverlapAlignment(
        TAlign & align,                         // [OUT] The target alignment object
        CandidateCoreSegmentMatch const & ccsm, // [IN]  The candidate core segment match
        TReadSequence & readSeq,
        unsigned const segId,                   // [IN]  The gene segment ID
        CdrReferences const & references,       // [IN]  The references
        double const & maxErrRate,              // [IN]  The maximum error rate allowed
        int shift,                              // [IN]  The SCF shift / offset
        LeftOverlap const                       // [TAG] The overlap direction
        )
{
    typedef typename Row<TAlign>::Type                          TRow;
    typedef typename Position<TRow>::Type                       TRowPos;
    typedef typename Source<TRow>::Type                         TSource;
    typedef typename Position<TSource>::Type                    TPos;
    typedef Segment<String<Dna5> const, InfixSegment>           TSegment;

    // Prepare the target align object
    resize(rows(align), 2);

    // The segment infix ends with the motif
    TPos segmentEndPos = references.leftMeta[segId].motifPos + 3;
    // If the SCF is shifted into the CDR3 region, extend
    if (shift > 0)
        segmentEndPos += shift;
    // Build the segments
    TSegment segSegment(references.leftSegs[segId], 0, segmentEndPos);
    TSegment readSegment(readSeq); // Dummy infix

    unsigned maxErrors = static_cast<unsigned>(std::ceil(1.0 * ccsm.readEndPos * maxErrRate));

    setSource(row(align, 0), readSegment);
    setSource(row(align, 1), segSegment);
    detach(align); // Copy infix but not data

    int diag = - length(segSegment) + ccsm.readEndPos;
    if (shift < 0)
        diag += -shift;
    if (diag > 0)
        diag = 0;

    // Compute the overlap alignment - we allow extra read sequence here
    int s = globalAlignment(align, SimpleScore(1,-1,-1), AlignConfig<true,true,false,true>(), diag - maxErrors, diag + maxErrors + 1);

    // Remove trailing gaps in the read sequence and pull in mismatches
    // instead. The score is equivalent, but we don't expect gaps here.
    refineTerminalGaps(align);

    // Clip after segment end
    TRowPos clippedViewEndPos = toViewPosition(row(align, 1), (length(segSegment)-1)) + 1;
    setClippedEndPosition(row(align, 0), clippedViewEndPos);
    setClippedEndPosition(row(align, 1), clippedViewEndPos);

    // Clip leading gaps
    size_t lGaps = std::max(countLeadingGaps(row(align, 0)), countLeadingGaps(row(align, 1)));
    setClippedBeginPosition(row(align, 0), lGaps);
    setClippedBeginPosition(row(align, 1), lGaps);

    return s;
}

template <typename TAlign, typename TReadSequence>
int extendToOverlapAlignment(
        TAlign & align,                         // [OUT] The target alignment object
        CandidateCoreSegmentMatch const & ccsm, // [IN]  The candidate core segment match
        TReadSequence & readSeq,
        unsigned const segId,                   // [IN]  The gene segment ID
        CdrReferences const & references,       // [IN]  The references
        double const & maxErrRate,              // [IN]  The maximum error rate allowed
        int shift,                              // [IN]  The SCF shift / offset
        RightOverlap const                      // [TAG] The overlap direction
        )
{
    typedef typename Row<TAlign>::Type                          TRow;
    typedef typename Position<TRow>::Type                       TRowPos;
    typedef typename Source<TRow>::Type                         TSource;
    typedef typename Position<TSource>::Type                    TPos;
    typedef Segment<String<Dna5> const, InfixSegment>           TSegment;

    // Prepare the target align object
    resize(rows(align), 2);

    // The segment infix begins with the motif
    TPos segmentBeginPos = references.rightMeta[segId].motifPos;
    // If the SCF is shifted into the CDR3 region, extend
    if (shift < 0)
        segmentBeginPos += shift;
    // Build the segments
    TSegment segSegment(references.rightSegs[segId], segmentBeginPos, length(references.rightSegs[segId]));
    TSegment readSegment(readSeq); // Dummy infix

    unsigned maxErrors = static_cast<unsigned>(length(readSeq) - ccsm.readBeginPos);

    setSource(row(align, 0), readSegment);
    setSource(row(align, 1), segSegment);
    detach(align); // Copy infix but not data

    int diag = ccsm.readBeginPos;

    // Compute the overlap alignment - we allow extra read sequence here
    int s = globalAlignment(align, SimpleScore(1,-1,-1), AlignConfig<true,false,true,true>(), diag - maxErrors, diag + maxErrors + 1);

    // Remove trailing gaps in the read sequence and pull in mismatches
    // instead. The score is equivalent, but we don't expect gaps here.
    refineTerminalGaps(align);

    // Clip after read or segment ends
    TRowPos readViewEndPos = toViewPosition(row(align, 0), length(readSegment) - 1) + 1;
    TRowPos segViewEndPos = toViewPosition(row(align, 1), length(segSegment) - 1) + 1;
    TRowPos clipViewEndPos = readViewEndPos > segViewEndPos ? segViewEndPos : readViewEndPos;
    setClippedEndPosition(row(align, 0), clipViewEndPos);
    setClippedEndPosition(row(align, 1), clipViewEndPos);

    // Clip before segment starts
    TRowPos clippedViewBeginPos = toViewPosition(row(align, 1), 0);
    setClippedBeginPosition(row(align, 0), clippedViewBeginPos);
    setClippedBeginPosition(row(align, 1), clippedViewBeginPos);

    return s;
}

inline double errRateFromScore(int score, unsigned length)
{
    return 1.0 * (length - score) / 2.0 / length;
}

/*
 * Identifies the best matching segment core fragment(s) among the provided
 * candidates
 */
template <typename TSegmentMatchesSet, typename TSequenceSet, typename TOverlapDirection>
void findBestSCFs(
        TSegmentMatchesSet & segMatchSet,                                       // [OUT] The resulting segment matches
        StringSet<String<CandidateCoreSegmentMatch> > const & candidateMatches, // [IN]  The candidate SCF-read matches
        TSequenceSet const & readSeqs,                                          // [IN]  The read sequences
        CdrReferences const & references,                                       // [IN]  The segment reference data
        CdrOptions const & options,                                             // [IN]  Runtime options
        std::vector<std::set<unsigned> > const * limSegmentIDs,                 // [IN]  Reduced sets of segment IDs to take into account
        TOverlapDirection const &                                               // [TAG] Indicating the overlap direction
        )
{
    // Type definitions
    typedef StringSet<String<CandidateCoreSegmentMatch> > const TCandMatchesSet;
    typedef Value<TCandMatchesSet>::Type const                  TCandMatches;
    typedef typename Value<TSequenceSet>::Type const            TSequence;
    typedef typename Infix<TSequence>::Type                     TInfix;
    typedef Align<TInfix>                                       TAlign;
    typedef typename Value<TSegmentMatchesSet>::Type            TSegmentMatches;
    typedef typename Value<TSegmentMatches>::Type               TSegmentMatch;

    TSequenceSet const & segSeqs = getSegmentSequences(references, TOverlapDirection());
    StringSet<String<unsigned> > const & scfToSegIds = getSCFToSegIds(references, TOverlapDirection());
    CdrReferences::TSCFPos const & scfBeginEndPos = getSCFPos(references, TOverlapDirection());

    // Clear output data
    clear(segMatchSet);
    resize(segMatchSet, length(readSeqs));

    double maxErrRate = getMaxErrRate(options, TOverlapDirection());

    // Iterate through the reads
    for (typename Iterator<TSequenceSet const, Rooted>::Type readIt = begin(readSeqs); !atEnd(readIt); goNext(readIt)) {
        TSequence  & readSeq    = *readIt;
        unsigned readId         = position(readIt);

        TSegmentMatches & segMatches = segMatchSet[readId];
        int maxScore = 0;
        // Iterate through the candidate matches and find the best one
        for (CandidateCoreSegmentMatch const & ccsm : candidateMatches[readId])
        {
            for (unsigned const segmentId : scfToSegIds[ccsm.coreSegId])
            {
                // Skip if we have a limiting set of ids and this one is not listed
                if (limSegmentIDs != nullptr && (*limSegmentIDs)[readId].find(segmentId) == (*limSegmentIDs)[readId].end())
                    continue;

                TAlign align;

                // Global overlap alignment computation
                int score = extendToOverlapAlignment(align, ccsm, readSeq, segmentId, references, maxErrRate, 
                        getSCFOffset(options, TOverlapDirection()), TOverlapDirection());

                double errRate = errRateFromScore(score, length(row(align, 0)));
                if (errRate > maxErrRate)
                    continue;

                // Store the result if there was no result before or if the
                // segment alignment is at least as good as the previous one(s)
                if (empty(segMatches) || maxScore <= score) {
                    if (!empty(segMatches) && maxScore < score) {
                        clear(segMatches);
                    }
                    maxScore = score;
                    appendValue(segMatches, TSegmentMatch(segmentId, score, align));
                }
            }
        }
    }
}

template <typename TSegmentMatchesSet, typename TSequenceSet, typename TOverlapDirection>
void findBestSCFs(
        TSegmentMatchesSet & segMatchSet,                                       // [OUT] The resulting segment matches
        StringSet<String<CandidateCoreSegmentMatch> > const & candidateMatches, // [IN]  The candidate SCF-read matches
        TSequenceSet const & readSeqs,                                          // [IN]  The read sequences
        CdrReferences const & references,                                       // [IN]  The segment reference data
        CdrOptions const & options,                                             // [IN]  Runtime options
        TOverlapDirection const &                                               // [TAG] Indicating the overlap direction
        )
{
    findBestSCFs(segMatchSet, candidateMatches, readSeqs, references, options, nullptr, TOverlapDirection());
}

/********************************************************************************
 * FUNCTION: findBestSegmentMatch()
 *
 * Wrapper function calling everything necessary to find the best segment matches
 * for either V or J segments
 *******************************************************************************/

template<typename TMatches, typename TGlobalData, typename TOverlapSpec>
void findBestSegmentMatch(
        TMatches & matches,                         // [OUT] Matches are stored here
        QueryData<SingleEnd> const & queryData,     //  [IN] The read sequences
        TGlobalData const & global,                 //  [IN] The segment references
        TOverlapSpec const)                         // [TAG] Overlap specification
{
    CdrReferences const & references = global.references;

    StringSet<String<CandidateCoreSegmentMatch> > candidateMatches;

    StringSet<TQueryDataSequence> const & seqs = getReadSequences(queryData, TOverlapSpec());

    CdrReferences::TSegCoreFragmentStringSet scfs = getSCFs(references, TOverlapSpec());

    findCandidateCoreSegments(candidateMatches, seqs, scfs, getMaxCoreSegErrors(global.options, TOverlapSpec()));

    findBestSCFs(
            matches,
            candidateMatches,
            seqs,
            references,
            global.options,
            TOverlapSpec());
}

template<typename TMatches, typename TGlobalData>
void findBestSegmentMatch(
        TMatches & matches,                         // [OUT] Matches are stored here
        QueryData<PairedEnd> const & queryData,     //  [IN] The read sequences
        TGlobalData const & global,                 //  [IN] The segment references
        RightOverlap const)                         // [TAG] Overlap specification
{
    CdrReferences const & references = global.references;

    StringSet<String<CandidateCoreSegmentMatch> > candidateMatches;

    StringSet<TQueryDataSequence>  const & seqs = getReadSequences(queryData, RightOverlap());

    CdrReferences::TSegCoreFragmentStringSet scfs = getSCFs(references, RightOverlap());

    findCandidateCoreSegments(candidateMatches, seqs, scfs, getMaxCoreSegErrors(global.options, RightOverlap()));

    findBestSCFs(
            matches,
            candidateMatches,
            seqs,
            references,
            global.options,
            RightOverlap());
}

template <typename TSequence, typename TGlobalData>
void findBestVSegment(
        std::vector<std::set<unsigned> > & dbMatches, // [OUT] For every read, the best matching V refs
        StringSet<TSequence> const & vReadSeqs,       //  [IN] V read sequences
        TGlobalData const & global)                   //  [IN] The segment references
{
    typedef typename Value<TSequence>::Type                             TAlphabet;

    // Types required for filtering
    typedef StringSet<TSequence>                                        TSequenceSet;
    typedef Shape<TAlphabet, SimpleShape>                               TShape;
    typedef Index<TSequenceSet, IndexQGram<TShape, OpenAddressing> >    TIndex;
    typedef Pattern<TIndex, Swift<SwiftSemiGlobal> >                    TPattern;
    typedef Finder<TSequence, Swift<SwiftSemiGlobal> >                  TFinder;

    // Types required for verification
    typedef Pattern<TSequence, MyersUkkonen>                            TMyersPattern;
    
    // Abort if no sequences were passed
    if (empty(vReadSeqs)) {
        dbMatches.clear();
        return;
    }

    // References (pattern)
    StringSet<TSequence> const & segmentSequences = global.references.leftSegs;

    double const maxErrRate = global.options.maxErrRateV;

    // ####################################################################################################
    // FILTERING
    // ####################################################################################################

    // Pattern is constructed over read sequences
    TIndex readsIndex(vReadSeqs);
    TPattern readsPattern(readsIndex);

    // Compute the q-gram length based on the desired error rate
    resize(indexShape(readsIndex), static_cast<int>(std::floor(1.0/maxErrRate)));

    // Here we store the potentially matching reads per segments
    String<std::set<unsigned> > readCandidates;
    resize(readCandidates, length(segmentSequences));

    for (typename Iterator<StringSet<TSequence> const, Rooted>::Type segSeqIt = begin(segmentSequences); !atEnd(segSeqIt); goNext(segSeqIt))
    {
        unsigned const segId = position(segSeqIt);
        TSequence & segSeq = const_cast<TSequence &>(*segSeqIt);
        TFinder segFinder(segSeq);

        while (find(segFinder, readsPattern, maxErrRate)) 
            readCandidates[segId].insert(segFinder.curHit->ndlSeqNo);
    }

    // ####################################################################################################
    // VERIFICATION
    // ####################################################################################################

    clear(dbMatches);
    dbMatches.resize(length(vReadSeqs));

    String<int> bestScores;
    resize(bestScores, length(vReadSeqs), MinValue<int>::VALUE);
    for (Iterator<String<std::set<unsigned> >, Rooted>::Type reCaIt = begin(readCandidates); !atEnd(reCaIt); goNext(reCaIt))
    {
        if (empty(*reCaIt))
            continue;

        unsigned const segId = position(reCaIt);
        
        for (std::set<unsigned>::const_iterator readIdIt = reCaIt->cbegin(); readIdIt != reCaIt->cend(); ++readIdIt)
        {

            unsigned vReadId = *readIdIt;

            TSequence const & vReadSeq = vReadSeqs[vReadId];
            int maxErrors = std::ceil(maxErrRate * length(vReadSeq));

            TMyersPattern myersPattern;
            setHost(myersPattern, vReadSeq);

            bool matchFound = false;
            int segBestScore = MinValue<int>::VALUE;

            Finder<TSequence const> finder(segmentSequences[segId]);
            while (find(finder, myersPattern, -maxErrors)) {
                if (!matchFound) {
                    segBestScore = getScore(myersPattern);
                    matchFound = true;
                } else {
                    int newScore = getScore(myersPattern);
                    if (newScore > segBestScore)
                        segBestScore = newScore;
                }
            }

            if (!matchFound) {
                continue;
            } 

            std::set<unsigned> & bestDBs = dbMatches[vReadId];

            if (segBestScore > bestScores[vReadId]) {
                bestDBs.clear();
                bestDBs.insert(segId);
                bestScores[vReadId] = segBestScore;
            } else if (segBestScore >= bestScores[vReadId]) {
                bestDBs.insert(segId);
            }
        }
    }
}

template<typename TMatches, typename TGlobalData>
void findBestSegmentMatch(
        TMatches & matches,                         // [OUT] Matches are stored here
        QueryData<PairedEnd> const & queryData,     //  [IN] The read sequences
        TGlobalData const & global,                 //  [IN] The segment references
        LeftOverlap const)                          // [TAG] Overlap specification
{
    typedef typename Value<TMatches>::Type              TSegmentMatches;
    typedef typename Value<TSegmentMatches>::Type       TSegmentMatch;
    typedef QueryData<PairedEnd>::TSequence             TSequence;
    typedef Segment<TSequence const, InfixSegment>      TSegment;
    typedef Align<TSegment>                             TAlign;

    typedef CdrReferences::TSegCoreFragmentStringSet                    TStringSet;
    typedef Value<TStringSet>::Type                                     TRefSequence;
    typedef typename Value<TRefSequence>::Type                          TAlphabet;
    typedef Shape<TAlphabet, SimpleShape>                               TShape;
    typedef Index<TStringSet, IndexQGram<TShape, OpenAddressing> >      TIndex;
    typedef Pattern<TIndex, Swift<SwiftSemiGlobal> >                    TPattern;
    typedef Finder<TSequence, Swift<SwiftSemiGlobal> >                  TFinder;
    typedef std::map<unsigned, BeginEndPos<long long> >                 TInfixPosMap;

    CdrReferences const & references = global.references;

    // Find V segments
    std::vector<std::set<unsigned> > vSegments;
    findBestVSegment(vSegments, queryData.fwSeqs, global);

    // Build index for SCF filtering
    unsigned maxSCFErrors = getMaxCoreSegErrors(global.options, LeftOverlap());
    TIndex scfIndex(references.leftSCFs);
    resize(indexShape(scfIndex), length(references.leftSCFs[0]) / (maxSCFErrors+1));
    TPattern segmentPattern(scfIndex);

    StringSet<TSequence> const & vSegSequences = references.leftSegs;

    // The length of all SCFs is the same
    double const errRate = 1.0 * global.options.maxVCoreErrors / length(references.leftSCFs[0]);

    // Iterate over all reads
    StringSet<String<CandidateCoreSegmentMatch> > scfResults;
    for (unsigned readId = 0; readId < vSegments.size(); ++readId)
    {
        TSequence & readSeq = const_cast<TSequence &>(queryData.revSeqs[readId]);

        std::set<unsigned> limScfIds;
        std::set<unsigned> const & limSegIds = vSegments[readId];

        if (limSegIds.empty())
        {
            appendValue(scfResults, String<CandidateCoreSegmentMatch>());
        }
        else
        {
            for (unsigned segId : limSegIds)
                limScfIds.insert(references.leftSegToScfId[segId]);

            // Filter
            TInfixPosMap scfInfixPositions = filterSCFs(readSeq, segmentPattern, errRate, &limScfIds);

            // Verification
            String<CandidateCoreSegmentMatch> ccsms = verifySCFHits(readSeq,
                    references.leftSCFs,
                    scfInfixPositions,
                    global.options.maxVCoreErrors);

            appendValue(scfResults, ccsms);
        }
    }

    // Find best segments
    StringSet<TQueryDataSequence> const & seqs = getVDJReadSequences(queryData);
    findBestSCFs(
            matches,
            scfResults,
            seqs,
            references,
            global.options,
            &vSegments,
            LeftOverlap());
}

#endif
