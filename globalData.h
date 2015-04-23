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

#include "sequence_data.h"

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
    typedef Shape<Dna5Q, SimpleShape> TShape;
    typedef StringSet<String<Dna5Q> >                                           TSegmentStringSet;
    typedef StringSet<String<Dna5Q> >                                           TSegCoreFragmentStringSet;
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
};

struct CdrOptions {
    CharString refFasta ;
    CharString rlogPath;
    CharString aminoOut;
    CharString nucOut;
    CharString fullOut;
    std::string outFileBaseName;
    unsigned vCrop;
    unsigned jCrop;
    int qmin;
    unsigned maxQClusterErrors;
    unsigned maxSClusterErrors;
    double maxClusterRatio;
    int maxQualityValue;
    int jobs;
    int pairedMinVOverlap;
    unsigned trunkReads;
    unsigned maxBlockSize;
    double minSdDevi;
    double qminclust;
    bool reverse;
    bool mergeAllels;
    bool cacheMatches;
    bool qualClustering;
    bool simpleClustering;
    bool mergeIdenticalCDRs;
    bool pairedEnd;
    bool bcRevRead;
    bool outputAligments;
    double maxErrRateV;
    double maxErrRateJ;
    double pairedMaxErrRateVOverlap;
    unsigned maxVCoreErrors;
    unsigned maxJCoreErrors;
    unsigned vSCFLength;
    unsigned jSCFLength;
    int vSCFOffset;
    int jSCFOffset;
    bool vSCFLengthAuto;
    unsigned vReadCrop;
    
    CdrOptions() : jobs(1), reverse(false), mergeAllels(false), cacheMatches(false), qualClustering(false), simpleClustering(false), mergeIdenticalCDRs(false), pairedEnd(false), bcRevRead(false), maxErrRateV(0), maxErrRateJ(0), maxVCoreErrors(0), maxJCoreErrors(0), vSCFLength(0), jSCFLength(0), vSCFOffset(-999), jSCFOffset(-999), vSCFLengthAuto(false), vReadCrop(0) {}
};

struct CdrOutputFiles {
    std::ofstream*      _fullOutStream;
    ConditionalLog      clusterCLog;
    Log                 clusterEvalLog;
};

template<typename TSequencingType>
struct CdrInputStreams {};

template<>
struct CdrInputStreams<SingleEnd> {
    std::string path;
    SeqFileIn stream;
    CdrInputStreams<SingleEnd>(std::string path_) : path(path_) {
        open(stream, path.c_str());
    }
};

template<>
struct CdrInputStreams<PairedEnd> {
    std::string fwPath, revPath;
    SeqFileIn fwStream, revStream;
    CdrInputStreams<PairedEnd>(std::string fwPath_, std::string revPath_) : fwPath(fwPath_), revPath(revPath_) {
        open(fwStream,fwPath.c_str());
        open(revStream,revPath.c_str());
    }
};

template<typename TSequencingType>
struct CdrGlobalData {
    CdrOptions const &                  options;
    CdrReferences const &               references;
    CdrInputStreams<TSequencingType> &  input;
    CdrOutputFiles &                    outFiles;

    CdrGlobalData(CdrOptions const & _options, CdrReferences const & _references, CdrInputStreams<TSequencingType> & _input, CdrOutputFiles & _outFiles) : 
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
// Getter for V/J segment alignment parameters
// ============================================================================

inline unsigned getMaxCoreSegErrors(CdrOptions const & options, LeftOverlap const)
{
    return options.maxVCoreErrors;
}

inline unsigned getMaxCoreSegErrors(CdrOptions const & options, RightOverlap const)
{
    return options.maxJCoreErrors;
}

inline double getMaxSegErrorRate(CdrOptions const & options, LeftOverlap const)
{
    return options.maxErrRateV;
}

inline double getMaxSegErrorRate(CdrOptions const & options, RightOverlap const)
{
    return options.maxErrRateJ;
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

// ============================================================================
// Getter for maximum error rate
// ============================================================================

inline double getMaxErrRate(CdrOptions const & options, LeftOverlap const)
{
    return options.maxErrRateV;
}

inline double getMaxErrRate(CdrOptions const & options, RightOverlap const)
{
    return options.maxErrRateJ;
}

#endif
