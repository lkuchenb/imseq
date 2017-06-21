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

#ifndef CDR3FINDER_GLOBAL_DATA_H
#define CDR3FINDER_GLOBAL_DATA_H

#include <iostream>

#include <seqan/index.h>
#include <seqan/arg_parse.h>

#include "runtime_options.h"
#include "sequence_data.h"
#include "segment_meta.h"
#include "logging.h"
#include "fastq_io.h"

using namespace seqan;

// ============================================================================
// CLASSES AND OPERATORS
// ============================================================================


struct ConditionalLog {
    std::ofstream* ofs;
    ConditionalLog() : ofs(NULL) {};
};

template<typename T>
ConditionalLog const & operator<<(ConditionalLog const & a, T const & b) {
    if (a.ofs!=NULL)
        (*a.ofs) << b;
    return a;
}

typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
typedef CoutType& (*ostreamFun)(CoutType&);

inline ConditionalLog const & operator<<(ConditionalLog const & a, ostreamFun fun)
{
    if (a.ofs!=NULL)
        fun(*(a.ofs));
    return a;
}

/**
 * Pair containing a begin and an end position
 */
template<typename TPos>
struct BeginEndPos {
    TPos beginPos, endPos;
    BeginEndPos() : beginPos(-1u), endPos(-1u) {}
    BeginEndPos(TPos _beginPos, TPos _endPos) : beginPos(_beginPos), endPos(_endPos) {}
};

struct CdrReferences {
    typedef Shape<Dna5, SimpleShape> TShape;
    typedef StringSet<String<Dna5> >                                           TSegmentStringSet;
    typedef StringSet<String<Dna5> >                                           TSegCoreFragmentStringSet;
    typedef StringSet<String<unsigned> >                                        TSCFToSegIds;
    typedef Index<TSegmentStringSet, IndexQGram<TShape, OpenAddressing> >       TQGramIndex;
    typedef String<SegmentMeta>                                                 TSegmentMetas;
    typedef String<BeginEndPos<unsigned> >                                      TSCFPos;

    TSegmentStringSet rightSegs, leftSegs;                      // The sequences of the segments
    TSegCoreFragmentStringSet rightSCFs, leftSCFs;              // The sequences of the segment core fragments
    TSCFToSegIds leftSCFToSegIds, rightSCFToSegIds;             // For each core fragment the ids of the corresponding segments
    TQGramIndex leftIndex, rightIndex;                          // The q-gram indices for the core fragments
    String<SegmentMeta> leftMeta, rightMeta;                    // The meta information for the segments
    TSCFPos leftSCFPos, rightSCFPos;                            // For each segment the begin and end position of the core fragment
    StringSet<String<unsigned> > leftIdentOffsets,              // The left and right identity offsets
        rightIdentOffsets;
    String<unsigned> leftToFirstAllel, rightToFirstAllel;       // Map pointing to the ID of the first allel for all allels
    String<unsigned> leftSegToScfId, rightSegToScfId;           // Map from segment ID to SCF id
};

struct CdrOutputFiles {
    std::ofstream*      _fullOutStream;
    ConditionalLog      clusterCLog;
    Log                 clusterEvalLog;
};

template<typename TSequencingType>
struct CdrGlobalData {
    CdrOptions const &                  options;
    CdrReferences const &               references;
    SeqInputStreams<TSequencingType> &  input;
    CdrOutputFiles &                    outFiles;

    CdrGlobalData(CdrOptions const & _options, CdrReferences const & _references, SeqInputStreams<TSequencingType> & _input, CdrOutputFiles & _outFiles) : 
        options(_options),
        references(_references),
        input(_input),
        outFiles(_outFiles)
    {}
};

template<typename T>
struct Sequencing {
};

template<typename T>
struct Sequencing<CdrGlobalData<T> > {
    typedef T Type;
};

// ============================================================================
// FUNCTIONS
// ============================================================================

inline void setConditionalLog(ArgumentParser & parser, ConditionalLog & condLog, std::string const & shortName) {
    if (isSet(parser, shortName)) {
        CharString logPath;
        getOptionValue(logPath, parser, shortName);
        condLog.ofs = new std::ofstream(toCString(logPath));
        if (!condLog.ofs->good()) {
            std::cerr << "Could not open the logfile path '" << logPath << "'!" << std::endl;
            exit(1);
        }
    }
}


// ============================================================================
// Getter for segment sequences
// ============================================================================

inline CdrReferences::TSegmentStringSet & getSegmentSequences(CdrReferences & references, RightOverlap const)
{
    return references.rightSegs;
}

inline CdrReferences::TSegmentStringSet & getSegmentSequences(CdrReferences & references, LeftOverlap const)
{
    return references.leftSegs;
}

inline CdrReferences::TSegmentStringSet const & getSegmentSequences(CdrReferences const & references, RightOverlap const)
{
    return references.rightSegs;
}

inline CdrReferences::TSegmentStringSet const & getSegmentSequences(CdrReferences const & references, LeftOverlap const)
{
    return references.leftSegs;
}

// ============================================================================
// Getter for scfToSegIds
// ============================================================================

inline CdrReferences::TSCFToSegIds & getSCFToSegIds(CdrReferences & references, RightOverlap const)
{
    return references.rightSCFToSegIds;
}

inline CdrReferences::TSCFToSegIds & getSCFToSegIds(CdrReferences & references, LeftOverlap const)
{
    return references.leftSCFToSegIds;
}

inline CdrReferences::TSCFToSegIds const & getSCFToSegIds(CdrReferences const & references, RightOverlap const)
{
    return references.rightSCFToSegIds;
}

inline CdrReferences::TSCFToSegIds const & getSCFToSegIds(CdrReferences const & references, LeftOverlap const)
{
    return references.leftSCFToSegIds;
}

// ============================================================================
// Getter for SCFPos
// ============================================================================

inline String<BeginEndPos<unsigned int> > & getSCFPos(CdrReferences & references, RightOverlap const)
{
    return references.rightSCFPos;
}

inline String<BeginEndPos<unsigned int> > & getSCFPos(CdrReferences & references, LeftOverlap const)
{
    return references.leftSCFPos;
}

inline String<BeginEndPos<unsigned int> > const & getSCFPos(CdrReferences const & references, RightOverlap const)
{
    return references.rightSCFPos;
}

inline String<BeginEndPos<unsigned int> > const & getSCFPos(CdrReferences const & references, LeftOverlap const)
{
    return references.leftSCFPos;
}

// ============================================================================
// Getter for SCFs
// ============================================================================

inline CdrReferences::TSegCoreFragmentStringSet & getSCFs(CdrReferences & references, RightOverlap const)
{
    return references.rightSCFs;
}

inline CdrReferences::TSegCoreFragmentStringSet & getSCFs(CdrReferences & references, LeftOverlap const)
{
    return references.leftSCFs;
}

inline CdrReferences::TSegCoreFragmentStringSet const & getSCFs(CdrReferences const & references, RightOverlap const)
{
    return references.rightSCFs;
}

inline CdrReferences::TSegCoreFragmentStringSet const & getSCFs(CdrReferences const & references, LeftOverlap const)
{
    return references.leftSCFs;
}


#endif
