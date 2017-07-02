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

#ifndef CDR3FINDER_REFERENCE_PREPARATION_H
#define CDR3FINDER_REFERENCE_PREPARATION_H

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include "segment_meta.h"
#include "globalData.h"

template<typename TSequence>
void computeIdentOffsets(
        StringSet<String<unsigned> > & identOffsets,    // [OUT] The ident offsets
        StringSet<TSequence> const & segSequences,      //  [IN] The segment sequences
        String<SegmentMeta> const & segMetaInfo         //  [IN] The segment meta information
        )
{
    typedef StringSet<TSequence>                                TSequenceSet;
    typedef typename Iterator<TSequence const, Rooted>::Type    TConstSeqIt;

    clear(identOffsets);
    resize(identOffsets, length(segSequences));
    for (unsigned i=0; i<length(identOffsets); ++i) 
        resize(identOffsets[i], length(segSequences));

    for (typename Iterator<TSequenceSet const, Rooted>::Type seqIt1 = begin(segSequences); !atEnd(seqIt1); goNext(seqIt1)) {
        typename Iterator<TSequenceSet const, Rooted>::Type seqIt2 = seqIt1;
        for (goNext(seqIt2); !atEnd(seqIt2); goNext(seqIt2)) {
            std::pair<unsigned, unsigned> idPair(position(seqIt1), position(seqIt2));
            TConstSeqIt charIt1 = iter(*seqIt1, segMetaInfo[idPair.first].motifPos);
            TConstSeqIt charIt2 = iter(*seqIt2, segMetaInfo[idPair.second].motifPos);
            unsigned offset = 0;
            for (;; ++offset) {
                if (*charIt1 != *charIt2)
                    break;
                if (charIt1 == begin(*seqIt1) || charIt2 == begin(*seqIt2))
                    break;
                --charIt1;
                --charIt2;
            }
            identOffsets[idPair.first][idPair.second] = offset;
            identOffsets[idPair.second][idPair.first] = offset;
        }
    }
}

/*
 * Parses a fasta file containing left or right segment reference sequences.
 * The fasta file must contain sequence IDs that conform to the format
 * Gene-Name|Segment-Type|Segment-ID|Allel-ID|Motif-Pos
 *
 * Example:
 * >TRB|J|1-6|02|22
 * CTCCTATAATTCACCCCTCCACTTTGGGAACGGGACCAGGCTCACTGTGACAG
 *
 * where 22 is the begin-position of the Phe118-encoding triplet.
 */
template <typename TSequence>
bool loadSegmentFiles(
        StringSet<TSequence>& vSegSequences,             // Output: The segment sequences
        String<SegmentMeta> & vSegMetaInfo,              // Output: The segment meta information
        StringSet<String<unsigned> > vIdentOffsets,      // Output: The identity offsets
        StringSet<TSequence>& jSegSequences,             // Output: The segment sequences
        String<SegmentMeta> & jSegMetaInfo,              // Output: The segment meta information
        StringSet<String<unsigned> > jIdentOffsets,      // Output: The identity offsets
        CharString const & path)                         // Input: Path to the fasta file
{
    typedef StringSet<CharString>               TIDSeqSet;

    // ============================================================================
    // Try to open the input file
    // ============================================================================

    SeqFileIn inFile;
    if (!open(inFile, toCString(path))) {
        std::cerr << "Error opening reference file " << path << std::endl;
        return false;
    }

    clear(vSegSequences);
    clear(vSegMetaInfo);
    clear(jSegSequences);
    clear(jSegMetaInfo);

    // ============================================================================
    // Try to parse the FASTA file
    // ============================================================================

    StringSet<TSequence> sequences;
    TIDSeqSet            ids;
    try {
        readRecords(ids, sequences, inFile);
    } catch (Exception const & e) {
        std::cerr << "ERROR parsing reference file '" << path << "': " << e.what() << std::endl;
        return false;
    }

    // ============================================================================
    // Parse the meta information and subtract the crop offset if necessary
    // ============================================================================

    for (Iterator<TIDSeqSet, Rooted>::Type it = begin(ids); !atEnd(it); goNext(it)) {
        SegmentMeta meta;
        if (!parseMetaInformation(meta, *it)) {
            std::cerr << "Invalid meta information (" << *it  << ") in file " << path << std::endl;
            return false;
        }
        if (meta.segType == "V") {
            appendValue(vSegSequences, sequences[position(it)]);
            appendValue(vSegMetaInfo, meta);
        } else if (meta.segType == "J") {
            appendValue(jSegSequences, sequences[position(it)]);
            appendValue(jSegMetaInfo, meta);
        } else if (meta.segType == "D") {
            continue;
        } else {
            std::cerr << "Invalid meta information (" << *it  << ") in file " << path << std::endl;
            return false;
        }
    }

    // ============================================================================
    // For each pair of sequences, identify the rightmost position left of the
    // motif where they differ and store it as an offset from the motif position
    // ============================================================================

    computeIdentOffsets(vIdentOffsets, vSegSequences, vSegMetaInfo);
    computeIdentOffsets(jIdentOffsets, jSegSequences, jSegMetaInfo);

    return true;
}

void buildToFirstAllelMap(String<SegmentMeta> const & meta, String<unsigned> & toFirstAllelMap)
{
    std::map<CharString, unsigned> idMap;
    clear(toFirstAllelMap);
    resize(toFirstAllelMap, length(meta));
    for (Iterator<String<SegmentMeta> const, Rooted>::Type it = begin(meta); !atEnd(it); goNext(it))
    {
        std::map<CharString, unsigned>::iterator oldIdIt = idMap.find(it->segId);
        if (oldIdIt == idMap.end())
        {
	    idMap[it->segId] = position(it);
            toFirstAllelMap[position(it)] = position(it);
        }
        else
        {
            toFirstAllelMap[position(it)] = oldIdIt->second;
        }
    }
}

void readAndPreprocessReferences(CdrReferences & references, CdrOptions const & options) {

    // ============================================================================
    // Read the segment file
    // ============================================================================

    if (!loadSegmentFiles(references.leftSegs, references.leftMeta, references.leftIdentOffsets, references.rightSegs, references.rightMeta, references.rightIdentOffsets, options.refFasta)) {
        std::cerr << "Reading the segments reference sequences failed!" << std::endl;
        exit(1);
    }

    // ============================================================================
    // PROCESS J-SEGMENTS
    // ============================================================================

    unsigned val = buildSegmentCoreFragments(references.rightSCFs, references.rightSCFToSegIds, references.rightSegToScfId, references.rightSCFPos, references.rightSegs, references.rightMeta, options.jSCFOffset, options.jSCFLength);
    if (BUILD_SEGMENT_CORE_FRAGMENTS_GOOD != val) {
        std::cerr << "Failed to build core fragment for '" << getDescriptor(references.rightMeta[val]) << "'. Boundaries violated.\n";
        exit(1);
    }

    buildToFirstAllelMap(references.rightMeta, references.rightToFirstAllel);

    // ============================================================================
    // PROCESS V-SEGMENTS
    // ============================================================================

    val = buildSegmentCoreFragments(references.leftSCFs, references.leftSCFToSegIds, references.leftSegToScfId, references.leftSCFPos, references.leftSegs, references.leftMeta, 3 - static_cast<int>(options.vSCFLength) + options.vSCFOffset, options.vSCFLength);
    if (BUILD_SEGMENT_CORE_FRAGMENTS_GOOD != val) {
        std::cerr << "Failed to build core fragment for '" << getDescriptor(references.leftMeta[val]) << "'. Boundaries violated.\n";
        exit(1);
    }

    buildToFirstAllelMap(references.leftMeta, references.leftToFirstAllel);

}
#endif
