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

// ############################################################################
// ############################################################################
// DEBUG SECTION

/* ENABLE TO GET VERBOSE OUTPUT */
//#define FULLDEBUG

/* ENABLE TO GET THE RESULTS FOR EACH READ PRINTED TO STDOUT */
//#define FULLREADOUT

// END DEBUG SECTION
// ############################################################################
// ############################################################################


#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>

#if SEQAN_HAS_ZLIB
#include <zlib.h>
#endif  // #if SEQAN_HAS_ZLIB

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

std::ofstream *__rejectLog;
#define REJECTLOG if (__rejectLog != NULL) *__rejectLog

std::ofstream *__phredStatStream;
#define PHREDSTAT if (__phredStatStream != NULL) *__phredStatStream

#define POINTERSTREAM(S) if (S != NULL) *S

#include "cdr3_cli.h"
#include "imseq.h"
#include "timeFormat.h"
#include "referencePreparation.h"

using namespace seqan;

bool checkSameSegment(String<SegmentMeta> const & meta, unsigned i, unsigned j)
{
    return (meta[i].segId == meta[j].segId);
}

// Program entry point
int main(int argc, char ** argv)
{
    // ============================================================================
    // Parse the command line options
    // ============================================================================

    CdrOptions options;
    CdrReferences references;
    CdrOutputFiles outFiles;
    String<std::string> inFilePaths;

    parseCommandLine(options, inFilePaths, argc, const_cast<char const **>(argv));

    // ============================================================================
    // Print the program call for potential logging
    // ============================================================================

    std::cerr << "===== Program call\n";
    for (int i=0; i<argc; ++i)
	std::cerr << argv[i] << ' ';
    std::cerr << std::endl;

    // ============================================================================
    // Call the workflow
    // ============================================================================

    InputInformation inputInformation;

    if (options.pairedEnd)
    {
        CdrInputStreams<PairedEnd> is(
                inFilePaths[0],
                inFilePaths[1]);
        CdrGlobalData<PairedEnd> global(
                options,
                references,
                is,
                outFiles
                );
        return main_generic(inputInformation, global, options, references);
    } else {
        CdrInputStreams<SingleEnd> is(
                inFilePaths[0]);
        CdrGlobalData<SingleEnd> global(
                options,
                references,
                is,
                outFiles
                );
        return main_generic(inputInformation, global, options, references);
    }

}
