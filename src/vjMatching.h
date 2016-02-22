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

typedef Segment<String<Dna5Q> const, InfixSegment>      _Dna5QInfixSegment;
struct CandidateCoreSegmentMatch : RefMatch<_Dna5QInfixSegment> {
    unsigned coreSegId;
    CandidateCoreSegmentMatch(unsigned coreSegId_, int score_, Align<_Dna5QInfixSegment> align_) : RefMatch<_Dna5QInfixSegment>(score_, align_), coreSegId(coreSegId_) { }
};

/********************************************************************************
 * FUNCTIONS
 *******************************************************************************/

/**
 * Constructs an ordered container with dependant infix segments. The infixes
 * have length 'length' and start at position 'start' relative to the motif
 * position
 */
const unsigned BUILD_SEGMENT_CORE_FRAGMENTS_GOOD = -1u;
template<typename TStringSet>
unsigned buildSegmentCoreFragments(
        TStringSet & scfSet,                                    // [OUT] Set of segment core fragments
        StringSet<String<unsigned> > & scfToSegId,              // [OUT] Set of segment ids for each core fragment
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
    reserve(scfSet, resMap.size());
    reserve(scfToSegId, resMap.size());
    for (typename TResMap::const_iterator resMapIt = resMap.begin(); resMapIt!=resMap.end(); ++resMapIt) {
        appendValue(scfSet, resMapIt->first);
        appendValue(scfToSegId, resMapIt->second);
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

template <typename TStringSet>
void findCandidateCoreSegments(
        StringSet<String<CandidateCoreSegmentMatch> > & results,        // [OUT] One String of candidate SCFs for each input read sequence
        TStringSet const & readSequences,                               // [IN]  The input read sequences
        TStringSet const & segmentSequences,                            // [IN]  The core fragments of the segment sequences
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
    typedef Finder<TSequence, Swift<SwiftSemiGlobal> >                  TFinder;
    typedef typename Infix<TSequence const>::Type                       TInfix;
    typedef Align<TInfix>                                               TAlign;
    typedef std::map<unsigned, BeginEndPos<long long> >                 TInfixPosMap;

    // Constants
    Score<int,Simple> const sc(1,-1,-1);

    // Pattern is constructed over segment core fragments, Finder is called on read sequences
    double const errRate = 1.0 * maxErrors / length(segmentSequences[0]); // The length of all SCFs is the same 
    TIndex segmentIndex(segmentSequences);
    resize(indexShape(segmentIndex), length(segmentSequences[0]) / (maxErrors+1));
    TPattern segmentPattern(segmentIndex);

    clear(results);
    reserve(results, length(readSequences));

    // Iterate through the read sequences and find potential core fragment matches
    for (typename Iterator<TStringSet const, Rooted>::Type readIt = begin(readSequences); !atEnd(readIt); goNext(readIt)) {
        // Need to give up const here. Maybe this won't be necessary some time in the 
        // distant future... 
        TSequence & readSeq = const_cast<TSequence &>(*readIt);
        TFinder readFinder(readSeq);

        String<CandidateCoreSegmentMatch> ccsms;
        TInfixPosMap scfInfixPositions;

        // ####################################################################################################
        // FILTERING
        // ####################################################################################################

        // Identify the windows within the read we have to investigate for potential
        // core fragment matches
        while (find(readFinder, segmentPattern, errRate, 0)) {
            unsigned coreSegId  = readFinder.curHit->ndlSeqNo;
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

        // ####################################################################################################
        // VERIFICATION
        // ####################################################################################################
    
        typedef BeginEndPos<long long>  TBeginEndPos;

        // Find all positions where the SCF matches within the constraints
        // specified by the user (# of errors)
        for (TInfixPosMap::const_iterator scfInfixPosition = scfInfixPositions.begin(); scfInfixPosition != scfInfixPositions.end(); ++scfInfixPosition) {
            unsigned coreSegId = scfInfixPosition->first;
            TBeginEndPos const & infixPos = scfInfixPosition->second;
            // Build the infix-SCF alignment
            TInfix readInfix(readSeq, infixPos.beginPos, infixPos.endPos);
            // "Pseudoinfix" to maintain type equality
            TInfix segInfix(segmentSequences[coreSegId]);

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

            // BUILD ALIGN OBJECTS FOR ALL CANDIDATE REGIONS 

            for (Iterator<String<TBeginEndPos>,Rooted>::Type it = begin(toExpandInfixes); !atEnd(it); goNext(it))
            {
                // Build align object
                TAlign align;
                resize(rows(align), 2);
                assignSource(row(align, 0), infix(readInfix, it->beginPos, it->endPos));
                assignSource(row(align, 1), segInfix);
                detach(align); // Copies the infix but not the sequence

                int score = globalAlignment(align, sc, AlignConfig<true,false,false,true>());

                CandidateCoreSegmentMatch ccsm(coreSegId, score, align);
                appendValue(ccsms, ccsm);

            }
        }

        appendValue(results, ccsms);
    }
}

template <typename TAlign, typename TSequence, typename TBeginEndPos>
int extendToOverlapAlignment(
        TAlign & align,                         // [OUT] The target alignment object
        CandidateCoreSegmentMatch const & ccsm, // [IN]  The candidate core segment match
        TSequence const & readSeq,              // [IN]  The read sequence
        TSequence const & segSeq,               // [IN]  The segment sequence
        TBeginEndPos const & segSeqSCFPos,      // [IN]  The SCF position inside the segment
        double const & maxErrRate,              // [IN]  The maximum error rate allowed
        LeftOverlap const                       // [TAG] The overlap direction
        )
{
    typedef typename Row<TAlign>::Type                          TRow;
    typedef typename Position<TRow>::Type                       TRowPos;
    typedef typename Source<TRow>::Type                         TSource;
    typedef typename Position<TSource>::Type                    TPos;
    typedef Segment<TSequence const, InfixSegment>              TSegment;
    
    // Prepare the target align object
    resize(rows(align), 2);

    // Clip begin and end of the SCF match if there are leading / trailing
    // mismatches or gaps
    TAlign scfAlignment = ccsm.align;
    TRow & scfAlignReadRow = row(scfAlignment, 0);
    TRow & scfAlignSegRow = row(scfAlignment, 1);
    clipOuterGaps(scfAlignment);

    // Determine the last position of the read and segment that are covered by
    // the provided SCF alignment
    TPos lastSCFReadPos = endPosition(scfAlignReadRow) + beginPosition(source(scfAlignReadRow));
    TPos lastSCFSegPos = endPosition(scfAlignSegRow) + beginPosition(source(scfAlignSegRow)) + segSeqSCFPos.beginPos;

    // Prepare the segments and set up the align object
    TSegment readSegment(readSeq, 0, lastSCFReadPos);
    unsigned maxErrors = static_cast<unsigned>(std::ceil(length(readSegment) * maxErrRate)) + 1;
    unsigned targetSegLength = length(readSegment)+maxErrors;
    TPos segSegmentBeginPos = lastSCFSegPos > targetSegLength ? lastSCFSegPos - targetSegLength : 0;
    TSegment segSegment(segSeq, segSegmentBeginPos, lastSCFSegPos);
    setSource(row(align, 0), readSegment);
    setSource(row(align, 1), segSegment);
    detach(align);

    // Compute the overlap alignment
    int s = globalAlignment(align, SimpleScore(1,-1,-1), AlignConfig<false,true,false,true>(), -(2*maxErrors), 0);

    // Clip leading segment gaps
    TRowPos clippedViewBeginPos = toViewPosition(row(align, 0), 0);
    setClippedBeginPosition(row(align, 0), clippedViewBeginPos);
    setClippedBeginPosition(row(align, 1), clippedViewBeginPos);

    return s;
}

template <typename TAlign, typename TSequence, typename TBeginEndPos>
int extendToOverlapAlignment(
        TAlign & align,                         // [OUT] The target alignment object
        CandidateCoreSegmentMatch const & ccsm, // [IN]  The candidate core segment match
        TSequence const & readSeq,              // [IN]  The read sequence
        TSequence const & segSeq,               // [IN]  The segment sequence
        TBeginEndPos const & segSeqSCFPos,      // [IN]  The SCF position inside the segment
        double const & maxErrRate,              // [IN]  The maximum error rate allowed
        RightOverlap const                      // [TAG] The overlap direction
        )
{
    typedef typename Row<TAlign>::Type                          TRow;
    typedef typename Position<TRow>::Type                       TRowPos;
    typedef typename Source<TRow>::Type                         TSource;
    typedef typename Position<TSource>::Type                    TPos;
    typedef Segment<TSequence const, InfixSegment>              TSegment;
    
    // Prepare the target align object
    resize(rows(align), 2);

    // Clip begin and end of the SCF match if there are leading / trailing
    // mismatches or gaps
    TAlign scfAlignment = ccsm.align;
    TRow & scfAlignReadRow = row(scfAlignment, 0);
    TRow & scfAlignSegRow = row(scfAlignment, 1);
    clipOuterGaps(scfAlignment);

    // Determine the first position of the read and segment that are covered by
    // the provided SCF alignment
    TPos firstSCFReadPos = beginPosition(scfAlignReadRow) + beginPosition(source(scfAlignReadRow));
    TPos firstSCFSegPos = beginPosition(scfAlignSegRow) + beginPosition(source(scfAlignSegRow)) + segSeqSCFPos.beginPos;

    // Prepare the segments and set up the align object
    TSegment readSegment(readSeq, firstSCFReadPos, length(readSeq));
    unsigned maxErrors = static_cast<unsigned>(std::ceil(length(readSegment) * maxErrRate)) + 1;
    TPos segSegmentEndPos = std::min(length(segSeq), length(readSegment) + maxErrors+ firstSCFSegPos);
    TSegment segSegment(segSeq, firstSCFSegPos, segSegmentEndPos);
    setSource(row(align, 0), readSegment);
    setSource(row(align, 1), segSegment);
    detach(align);

    // Compute the overlap alignment
    int s = globalAlignment(align, SimpleScore(1,-1,-1), AlignConfig<true,false,true,true>(), -maxErrors, maxErrors);

    // Clip leading segment gaps
    TRowPos clippedViewEndPos = std::min(toViewPosition(row(align, 0), length(source(row(align,0)))), toViewPosition(row(align, 1), length(source(row(align,1)))));
    setClippedEndPosition(row(align, 0), clippedViewEndPos);
    setClippedEndPosition(row(align, 1), clippedViewEndPos);

    return s;
}

template <typename TAlign, typename TSequence, typename TBeginEndPos, typename TOverlapDirection>
int expandSCFAlignment(
        TAlign & align,                         // [OUT] The target alignment object
        CandidateCoreSegmentMatch const & ccsm, // [IN]  The candidate core segment match
        TSequence const & readSeq,              // [IN]  The read sequence
        TSequence const & segSeq,               // [IN]  The segment sequence
        TBeginEndPos const & segSeqSCFPos,      // [IN]  The SCF position inside the segment
        TOverlapDirection const &               // [TAG] The overlap direction
        )
{
    // Type definitions
    typedef typename Infix<TSequence const>::Type       TInfix;

    // The align object in the CandidateCoreSegmentMatch contains an Infix<->Infix alignment
    // between the segment core fragment (directly, not as infix of the segment) and the 
    // read. At first, the segment core fragment infix needs to be replaced by a segment infix.

    TInfix segInfix(segSeq, segSeqSCFPos.beginPos, segSeqSCFPos.endPos);
    align = ccsm.align;

    // Cannot use "setSource()" here because this would reset the gaps.
    // Maybe we should use a new alignment + integrateAlign() here instead?
    setValue(row(align, 1)._source, segInfix);

    detach(align);

    Tuple<unsigned, 4> positions;
    clipAndComputePositions(positions, align, TOverlapDirection());

    Score<int,Simple> const sc(1,-1,-1);

    int score = extendAlignment(align, ccsm.score, readSeq, segSeq, positions, GetExtDir<TOverlapDirection>::VALUE, -6, 6, sc);

    return score;
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
        double const maxErrRate,                                                // [IN]  The maximum error rate
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

    // Iterate through the reads
    for (typename Iterator<TSequenceSet const, Rooted>::Type readIt = begin(readSeqs); !atEnd(readIt); goNext(readIt)) {
        TSequence  & readSeq    = *readIt;
        unsigned readId         = position(readIt);

        TSegmentMatches & segMatches = segMatchSet[readId];
        int maxScore = 0;
        // Iterate through the candidate matches and find the best one
        for (typename Iterator<TCandMatches, Rooted>::Type cmIt = begin(candidateMatches[readId]); !atEnd(cmIt); goNext(cmIt)) {
            for (Iterator<String<unsigned> const, Rooted>::Type segIdIt = begin(scfToSegIds[cmIt->coreSegId]); !atEnd(segIdIt); goNext(segIdIt)) {
                unsigned segmentId = *segIdIt;

                TAlign align;

                // Global overlap alignment computation
                int score = extendToOverlapAlignment(align, *cmIt, readSeq, segSeqs[segmentId], scfBeginEndPos[segmentId], maxErrRate, TOverlapDirection());

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
            getMaxErrRate(global.options, TOverlapSpec()),
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
            getMaxErrRate(global.options, RightOverlap()),
            RightOverlap());
}

template <typename TSequence, typename TGlobalData>
void findBestVSegment(
        StringSet<String<unsigned> > & dbMatches, // [OUT] For every read, the best matching V refs
        StringSet<TSequence> const & vReadSeqs,   //  [IN] V read sequences
        TGlobalData const & global)               //  [IN] The segment references
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
        clear(dbMatches);
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
    resize(dbMatches, length(vReadSeqs));

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

            String<unsigned> & bestDBs = dbMatches[vReadId];

            if (segBestScore > bestScores[vReadId]) {
                clear(bestDBs);
                appendValue(bestDBs, segId);
                bestScores[vReadId] = segBestScore;
            } else if (segBestScore >= bestScores[vReadId]) {
                appendValue(bestDBs, segId);
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

    CdrReferences const & references = global.references;

    clear(matches);
    reserve(matches, nRecords(queryData));

    // (1) Find V segments

    StringSet<String<unsigned> > vSegments;
    findBestVSegment(vSegments, queryData.fwSeqs, global);

    // (3) Align V segment to VDJ read (overlap alignment)
    // TODO is an overlap alignment really necessary? Can we just make a semi-global
    // alignment, assuming that the read does not contain any non-V-segment sequence?

    StringSet<TSequence> const & vSegSequences = references.leftSegs;

    // Iterate over all reads
    for (Iterator<StringSet<String<unsigned> > const, Rooted>::Type it = begin(vSegments); !atEnd(it); goNext(it))
    {
        unsigned readId = position(it);
        TSegmentMatches readMatches;

        // Iterate over all matches
        int maxScore = INT_MIN;
        for (Iterator<String<unsigned> const, Rooted>::Type refIt = begin(*it); !atEnd(refIt); goNext(refIt))
        {
            // Compute the overlap alignment between the VDJ read and the V sequence
            TAlign align;
            resize(rows(align), 2);
            assignSource(row(align, 0), TSegment(queryData.revSeqs[readId]));

            TSegment segSegment(vSegSequences[*refIt]);

            // Clip off sequence that is downstream of the C104 triplet - this
            // limits the analysis and the alignment score / errRate to the
            // non-CDR3 region, just as in the single-end expansion case
            setEndPosition(segSegment, global.references.leftMeta[*refIt].motifPos+3);

            // Clip off extraneous sequence from the V-segment, that couldn't
            // be matched to to the limited read-length
            if (length(vSegSequences[*refIt]) > length(segSegment))
                setBeginPosition(segSegment, length(segSegment) - length(queryData.revSeqs[readId]));

            assignSource(row(align, 1), segSegment);

            // TODO Check hard-coded parameters
            int score           = globalAlignment(align, SimpleScore(1, -2, -2), AlignConfig<false,true,false,true>(), -length(segSegment), 10);

            // Clip to non-CDR3 overlap region
            int olBeginViewPos = toViewPosition(row(align, 0),0);
            int olEndViewPos = 1 + toViewPosition(row(align, 1), length(source(row(align,1))) - 1 );
            setClippedBeginPosition(row(align, 0), olBeginViewPos);
            setClippedBeginPosition(row(align, 1), olBeginViewPos);
            setClippedEndPosition(row(align, 0), olEndViewPos);
            setClippedEndPosition(row(align, 1), olEndViewPos);

            int overlapLength   = olEndViewPos - olBeginViewPos;
            int mismatches      = ((overlapLength - score) / 3);
            double errRate      = 1.0 * mismatches / overlapLength;
            int unitScore       = overlapLength - 2 * mismatches;

            // We keep track of the maximum observed score since we can still
            // discard alignments based on the V-segment against VDJ-read
            // alignment to refine the results obtained from the V-read
            // alignments
            if (unitScore < maxScore)
                continue;

            if (overlapLength >= global.options.pairedMinVOverlap && errRate <= global.options.pairedMaxErrRateVOverlap) {
                if (unitScore > maxScore) {
                    clear(readMatches);
                    maxScore = unitScore;
                }
                appendValue(readMatches, TSegmentMatch(*refIt, unitScore, align));
            }
        }

        appendValue(matches, readMatches);
    }
}

#endif
