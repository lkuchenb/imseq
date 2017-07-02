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

#ifndef CDRFINDER_CMDLINE_H
#define CDRFINDER_CMDLINE_H

#include <thread>

#include <seqan/arg_parse.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include "version_number.h"
#include "thread_check.h"
#include "logging.h"
#include "segment_meta.h"
#include "clone.h"
#include "vjMatching.h"
#include "globalData.h"

#define  OPT_QMIN_DEFAULT       30
#define  OPT_CQMIN_DEFAULT      0
#define  OPT_COPYJ_DEFAULT      0
#define  OPT_COPYV_DEFAULT      0
#define  OPT_EPSJ_DEFAULT       0.1
#define  OPT_EPSV_DEFAULT       0.1
#define  OPT_LEFT_MOTLEN        3u
#define  OPT_LEFT_MOTOFF        0
#define  OPT_LEFT_MOTERR        0u
#define  OPT_RIGHT_MOTLEN       3u
#define  OPT_RIGHT_MOTOFF       0
#define  OPT_RIGHT_MOTERR       0u
#define  OPT_TRUNC_READS        0u
#define  OPT_MAX_ERRORS_LQ      4
#define  OPT_MAX_ERRORS_HQ      2
#define  OPT_MAX_CLUST_RATIO    1.0
#define  OPT_CLUST_SD           1.0
#define  OPT_QMINCLUST_DEFAULT  30.0
#define  OPT_MAX_QUALITY_VALUE  40
#define  OPT_JOBS_DEFAULT       1
#define  OPT_VSCOREBOUNDARY     0
#define  OPT_JSCOREBOUNDARY     0
#define  OPT_VL_DEFAULT         10u
#define  OPT_JL_DEFAULT         10u
#define  OPT_MBL_DEFAULT        0

using namespace seqan;

inline void initializeLog(Log & log, ArgumentParser& parser, std::string const & shortName, std::string altName = std::string()) {
    std::string s;
    getOptionValue(s, parser, shortName);
    if (s=="-" && altName!="") s=altName;
    log.setPath(s);
}

inline void parseCommandLine(CdrOptions & options, String<std::string> & inFilePaths, const int argc, const char** argv) {

    ArgumentParser parser("imseq");
    setVersion(parser, IMSEQ_VERSION::STRING);
    setDate(parser, "March 2016");

    addUsageLine(parser, "-ref <segment reference> [\\fIOPTIONS\\fP] <VDJ reads>");
    addUsageLine(parser, "-ref <segment reference> [\\fIOPTIONS\\fP] <V reads> <VDJ reads>");

    addDescription(parser,
            "\\fBimseq\\fP is a tool for the analysis of T- and B-cell receptor chain sequences. It can be used "
            "to analyse either single-read data, where the reads cover the V-CDR3-J region "
            "sufficiently for an identification, or paired-end data where one read covers the "
            "V-region and one read covers the J- and CDR3-region. The latter read has do cover only a "
            "small fraction of the V-segment, sufficient for the localization of the Cys-104 motif.");

    addDescription(parser, "The following options exist:");

    // ============================================================================
    // One argument has to be specified, namely the input file
    // ============================================================================

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "<input file> [<input file>]", 1));

    // ============================================================================
    // Plenty of options
    // ============================================================================

    addSection(parser, "The following switch must be specified");
    addOption(parser, ArgParseOption("ref", "reference", "FASTA file with gene segment reference sequences.", ArgParseArgument::INPUT_FILE ));
    setRequired(parser, "ref");
    addSection(parser, "Output files. At least one of the following switches must be specified");
    addOption(parser, ArgParseOption("oa", "out-amino", "Output file path for translated clonotypes.", (ArgParseArgument::STRING)));
    addOption(parser, ArgParseOption("on", "out-nuc", "Output file path for untranslated clonotypes.", (ArgParseArgument::STRING)));
    addOption(parser, ArgParseOption("o", "out", "Output file path for verbose output per analyzed read.", (ArgParseArgument::STRING)));
    addOption(parser, ArgParseOption("s", "seq", "Include read sequence in output (-o) file."));

    //================================================================================
    // Read preprocessing
    //================================================================================

    addSection(parser, "Read preprocessing");
    addOption(parser, ArgParseOption("r", "reverse", "By default, V-reads are read as they are and V(D)J-reads are reverse complemented. Use this switch for the opposite behaviour."));
    addOption(parser, ArgParseOption("tr", "truncate-reads", "Truncate reads to the specified length. 0 leaves them at their original lengths.", (ArgParseArgument::INTEGER)));
    setDefaultValue(parser, "tr", OPT_TRUNC_READS);

    //================================================================================
    // Additional output
    //================================================================================

    addSection(parser, "Additional output settings");
    addOption(parser, ArgParseOption("rlog", "reject-log", "Log file for rejected reads. If empty, no log file is written.", (ArgParseArgument::OUTPUT_FILE)));
    addOption(parser, ArgParseOption("al", "with-alleles", "Keep allele information in output files and during aggregation."));

    //================================================================================
    // V / J segment alignment
    //================================================================================

    addSection(parser, "V/J segment alignment");
    addOption(parser, ArgParseOption("ev", "v-err-rate", "Maximum error rate allowed within the V segment alignment", ArgParseArgument::DOUBLE));
    setMinValue(parser, "ev", "0");
    setDefaultValue(parser, "ev", 0.05);
    addOption(parser, ArgParseOption("ej", "j-err-rate", "Maximum error rate allowed within the J segment alignment", ArgParseArgument::DOUBLE));
    setMinValue(parser, "ej", "0");
    setDefaultValue(parser, "ej", 0.15);

    //================================================================================
    // V / J segment alignment (paired end)
    //================================================================================

    addSection(parser, "V/J segment alignment (paired-end)");
    addOption(parser, ArgParseOption("pve", "paired-v-error", "Maximum error rate in the alignment between the forward-read identified V segment and the reverse read. Default: Use value from -ev.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "pve", "0");
    setMaxValue(parser, "pve", "1");
    addOption(parser, ArgParseOption("vcr", "v-read-crop", "Crop NUM bases from the beginning of the V read before processing it", ArgParseArgument::INTEGER));
    setMinValue(parser, "vcr", "0");
    setDefaultValue(parser, "vcr", 0);

    //================================================================================
    // V / J segment alignment (expert)
    //================================================================================

    addSection(parser, "V/J segment alignment (Expert settings)");
    addOption(parser, ArgParseOption("jcl", "j-core-length", "Length of the J core fragment.", ArgParseArgument::INTEGER));
    setMinValue(parser, "jcl", "5");
    setMaxValue(parser, "jcl", "20");
    setDefaultValue(parser, "jcl", 12);
    addOption(parser, ArgParseOption("jco", "j-core-offset", "Offset of the V core fragment.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "jco", -6);
    addOption(parser, ArgParseOption("vcl", "v-core-length", "Length of the V core fragment. Default: Automatically select value between 10 and 20 based on minimum observed read length.", ArgParseArgument::INTEGER));
    setMinValue(parser, "vcl", "5");
    addOption(parser, ArgParseOption("vco", "v-core-offset", "Offset of the V core fragment.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "vco", 0);
    addOption(parser, ArgParseOption("vce", "v-core-errors", "Maximum number of errors when matching the V core fragments.", ArgParseArgument::INTEGER));
    setMinValue(parser, "vce", "0");
    setDefaultValue(parser, "vce", 1);
    addOption(parser, ArgParseOption("jce", "j-core-errors", "Maximum number of errors when matching the J core fragments.", ArgParseArgument::INTEGER));
    setMinValue(parser, "jce", "0");
    setDefaultValue(parser, "jce", 2);

    //================================================================================
    // QUALITY CONTROL [Alignment score thresholds, q-score thresholds]
    //================================================================================

    addSection(parser, "Quality control");
    addOption(parser, ArgParseOption("mq", "min-qual", "Minimum average read phred score. In paired end mode, this is applied to both reads. See '-sfb'.", (ArgParseArgument::INTEGER)));
    setMinValue(parser, "mq", "0");
    setMaxValue(parser, "mq", "60");
    setDefaultValue(parser, "mq", OPT_QMIN_DEFAULT);
    addOption(parser, ArgParseOption("mcq", "min-clust-qual", "Minimum average cluster phred score.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "mcq", "0");
    setMaxValue(parser, "mcq", "60");
    setDefaultValue(parser, "mcq", OPT_QMIN_DEFAULT);
    addOption(parser, ArgParseOption("mrl", "min-read-length", "Minimum read length. In paired end mode, this is applied to both reads. See '-sfb'.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mrl", "0");
    setDefaultValue(parser, "mrl", 75);
    addOption(parser, ArgParseOption("mcl", "min-cdr3-length", "Minimum CDR3 length in amino acids.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mcl", "0");
    setDefaultValue(parser, "mcl", 5);
    addOption(parser, ArgParseOption("sfb", "single-end-fallback", "Fall back to single end analysis based on VDJ read if V read fails -mq or -mrl."));

    //================================================================================
    // BARCODING
    //================================================================================

    addSection(parser, "Barcoding");

    // SWITCHES
    addOption(parser, ArgParseOption("bvdj", "barcode-vdj", "In paired end mode: Read the barcode from the VDJ read instead of the V read."));
    // PARAMETERS
    addOption(parser, ArgParseOption("bse", "bcseq-max-err", "Maximum number of errors allowed in the barcode sequence", (ArgParseArgument::INTEGER)));
    setMinValue(parser, "bse", "0");
    setDefaultValue(parser, "bse", "1");
    addOption(parser, ArgParseOption("bmq", "bc-min-qual", "Minimum per base quality in molecular barcode region", (ArgParseArgument::INTEGER)));
    setMinValue(parser, "bmq", "0");
    setMaxValue(parser, "bmq", "60");
    setDefaultValue(parser, "bmq", "30");
    addOption(parser, ArgParseOption("bcl", "barcode-length", "Length of random barcode at the beginning of the read. A value of '0' disables barcode based correction.", (ArgParseArgument::INTEGER)));
    setMinValue(parser, "bcl", "0");
    setDefaultValue(parser, "bcl", 0);
    addOption(parser, ArgParseOption("ber", "barcode-err-rate", "Maximum error rate between reads in order to be merged based on barcode sequence", (ArgParseArgument::DOUBLE)));
    setMinValue(parser, "ber", "0");
    setMaxValue(parser, "ber", "1");
    setDefaultValue(parser, "ber", 0.05);
    addOption(parser, ArgParseOption("bfr", "barcode-freq-rate", "Inclusive maximum frequency ratio between smaller and larger cluster during barcode clustering", (ArgParseArgument::DOUBLE)));
    setMinValue(parser, "bfr", "0");
    setMaxValue(parser, "bfr", "1");
    setDefaultValue(parser, "bfr", 0.2);
    addOption(parser, ArgParseOption("bst", "barcode-stats", "Path to barcode stats output file. If empty, no file is written. Ignored if -bcl is 0.", (ArgParseArgument::OUTPUT_FILE)));
    addOption(parser, ArgParseOption("oab", "out-amino-bc", "Output file path for translated clonotypes with barcode corrected counts. Ignored if -bcl is 0.", (ArgParseArgument::OUTPUT_FILE)));
    addOption(parser, ArgParseOption("onb", "out-nuc-bc", "Output file path for untranslated clonotypes with barcode corrected counts. Ignored if -bcl is 0.", (ArgParseArgument::OUTPUT_FILE)));


    //================================================================================
    // POSTPROCESSING / CLUSTERING
    //================================================================================

    addSection(parser, "Postprocessing / Clustering");
    // SWITCHES
    addOption(parser, ArgParseOption("ma", "merge-ambiguous-seg", "Merge clonotypes with identical CDR3 sequences separated by ambiguous segment identification"));
    addOption(parser, ArgParseOption("qc", "qual-clust", "Enable quality score based clustering."));
    addOption(parser, ArgParseOption("sc", "simple-clust", "Enable simple distance-based clustering"));
    // PARAMETERS
    addOption(parser, ArgParseOption("qcme", "max-err-hq", "Maximum edit-distance for two clusters to be clustered without low quality correlation", (ArgParseArgument::INTEGER)));
    setMinValue(parser, "qcme", "0");
    setDefaultValue(parser, "qcme", OPT_MAX_ERRORS_LQ);
    addOption(parser, ArgParseOption("qcsd", "min-sd-diff", "How many standard deviations must an error positions quality value be below the mean to be considered for clustering.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "qcsd", OPT_CLUST_SD);
    addOption(parser, ArgParseOption("scme", "max-err-lq", "Maximum edit-distance for two clusters to be potentially clustered without low quality correlation", (ArgParseArgument::INTEGER)));
    setMinValue(parser, "scme", "0");
    setDefaultValue(parser, "scme", OPT_MAX_ERRORS_HQ);
    addOption(parser, ArgParseOption("mcr", "max-clust-ratio", "Maximum abundance ratio for two clonotypes to be clustered", (ArgParseArgument::DOUBLE)));
    setMinValue(parser, "mcr", "0");
    setMaxValue(parser, "mcr", "1");
    setDefaultValue(parser, "mcr", OPT_MAX_CLUST_RATIO);

    //================================================================================
    // PERFORMANCE [Parallelization, Caching]
    //================================================================================

    addSection(parser, "Performance");
#ifdef __WITHCDR3THREADS__
    addOption(parser, ArgParseOption("j", "jobs", "Number of parallel jobs (threads).", (ArgParseArgument::INTEGER)));
    setDefaultValue(parser, "j", OPT_JOBS_DEFAULT);
#endif

    //================================================================================
    // Other options
    //================================================================================


    addSection(parser, "Other options");
    addOption(parser, ArgParseOption("pa", "print-alignments", "Print the V/J alignments for each read. Implies -j 1."));

    // ============================================================================
    // Exit if there was an error or if the help switch was used
    // ============================================================================

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        exit(res == ArgumentParser::PARSE_ERROR);
    if (getArgumentValueCount(parser, 0) < 1 || getArgumentValueCount(parser, 0) > 2) {
        std::cerr << "You must specify exactly one or two input files!" << std::endl;
        exit(1);
    }

    // ============================================================================
    // Check some conditions that cannot be specified with the ArgumentParser
    // ============================================================================

    if (!isSet(parser, "oa") && !isSet(parser, "on") && !isSet(parser, "o")) {
        std::cerr << "You have to specify at least one of the following options: -o, -oa, -on\n";
        exit(1);
    }

    // ============================================================================
    // Read the command line argument determining the input filename
    // ============================================================================

    inFilePaths = getArgumentValues(parser, 0);

    options.pairedEnd = length(inFilePaths) == 2;

    SEQAN_CHECK(length(inFilePaths)==1 || length(inFilePaths)==2, "Please report this error");

    // ============================================================================
    // Determine the output basename
    // ============================================================================
    // -1- remove gzip suffix if present
    CharString outFile = inFilePaths[0];
    options.outFileBaseName = toCString(outFile);
    std::string testString = options.outFileBaseName.substr(options.outFileBaseName.size()-3);
    if (testString==".gz" || testString==".GZ" || testString==".Gz" || testString == ".gZ")
        options.outFileBaseName = options.outFileBaseName.substr(0, options.outFileBaseName.size()-3);
    // -2- remove the last '.' separated column
    size_t periodPos = options.outFileBaseName.rfind(".");
    if (periodPos != std::string::npos)
        options.outFileBaseName = options.outFileBaseName.substr(0, periodPos);

    // ============================================================================
    // Fill the options struct
    // ============================================================================

    getOptionValue(options.refFasta, parser, "ref");
    getOptionValue(options.aminoOut, parser, "oa");
    if (options.aminoOut=="-") options.aminoOut = options.outFileBaseName + ".act";
    getOptionValue(options.nucOut, parser, "on");
    if (options.nucOut=="-") options.nucOut = options.outFileBaseName + ".nct";
    getOptionValue(options.fullOut, parser, "o");
    if (options.fullOut=="-") options.fullOut = options.outFileBaseName + ".rdt";
    getOptionValue(options.rlogPath, parser, "rlog");
    if (options.rlogPath=="-") options.rlogPath = options.outFileBaseName + ".rlg";

    getOptionValue(options.barcodeLength, parser, "bcl");
    getOptionValue(options.bcClustMaxErrRate, parser, "ber");
    getOptionValue(options.bcClustMaxFreqRate, parser, "bfr");
    options.barcodeVDJRead = isSet(parser, "bvdj");
    options.rdtWithSequence = isSet(parser, "s");

    // -bst --barcode-stats
    if (isSet(parser, "bst") && options.barcodeLength > 0)
        getOptionValue(options.bstPath, parser, "bst");
    else
        options.bstPath = "";
    if (options.bstPath=="-") options.bstPath = options.outFileBaseName + ".bst";

    // -oab --out-amino-bc
    if (isSet(parser, "oab") && options.barcodeLength > 0)
        getOptionValue(options.aminoOutBc, parser, "oab");
    else
        options.aminoOutBc = "";
    if (options.aminoOutBc=="-") options.aminoOutBc = options.outFileBaseName + ".bact";

    // -onb --out-nuc-bc
    if (isSet(parser, "onb") && options.barcodeLength > 0)
        getOptionValue(options.nucOutBc, parser, "onb");
    else
        options.nucOutBc = "";
    if (options.nucOutBc=="-") options.nucOutBc = options.outFileBaseName + ".bnct";

    options.mergeAllels = !isSet(parser, "al"); // [!] Mind the negation
    options.reverse = isSet(parser, "r");
    getOptionValue(options.trunkReads, parser, "tr");
    options.maxQualityValue = 0;
    getOptionValue(options.maxErrRateV, parser, "ev");
    getOptionValue(options.maxErrRateJ, parser, "ej");
    getOptionValue(options.maxVCoreErrors, parser, "vce");
    getOptionValue(options.maxJCoreErrors, parser, "jce");

    getOptionValue(options.jSCFLength, parser, "jcl");
    getOptionValue(options.jSCFOffset, parser, "jco");
    getOptionValue(options.vSCFOffset, parser, "vco");
    options.vSCFLengthAuto = !isSet(parser, "vcl");
    if (!options.vSCFLengthAuto)
        getOptionValue(options.vSCFLength, parser, "vcl");

    getOptionValue(options.qmin, parser, "mq");
    getOptionValue(options.barcodeMaxError, parser, "bse");
    getOptionValue(options.bcQmin, parser, "bmq");
    getOptionValue(options.qminclust, parser, "mcq");
    getOptionValue(options.minReadLength, parser, "mrl");
    getOptionValue(options.minCDR3Length, parser, "mcl");
    if (isSet(parser, "pve"))
        getOptionValue(options.pairedMaxErrRateVOverlap, parser, "pve");
    else
        options.pairedMaxErrRateVOverlap = options.maxErrRateV;
    getOptionValue(options.vReadCrop, parser, "vcr");
    options.singleEndFallback = isSet(parser, "sfb");
    options.pairedMinVOverlap = 10;
    options.mergeIdenticalCDRs = isSet(parser, "ma");

    // CLUSTERING
    getOptionValue(options.maxQClusterErrors, parser, "qcme");
    getOptionValue(options.maxSClusterErrors, parser, "scme");
    getOptionValue(options.maxClusterRatio, parser, "mcr");
    options.qualClustering = isSet(parser, "qc") && (options.maxQClusterErrors > 0) ;
    options.simpleClustering = isSet(parser, "sc") && (options.maxSClusterErrors > 0);
    getOptionValue(options.minSdDevi, parser, "qcsd");
//    initializeLog(outFiles.clusterEvalLog, parser, "cvo", options.outFileBaseName + ".cel");

    options.cacheMatches = true;
#ifdef __WITHCDR3THREADS__
    if (isSet(parser, "j"))
    {
        getOptionValue(options.jobs, parser, "j");
    }
    else
    {
        options.jobs = std::thread::hardware_concurrency();
    }
#endif
//    setConditionalLog(parser, outFiles.clusterCLog, "cl");
    options.outputAligments = isSet(parser, "pa");
    if (options.outputAligments)
        options.jobs = 1;

    options.vCrop = options.pairedEnd ? 0 : 150;
    options.jCrop = 150;

    //TODO Does it make sense to make this configurable?
    options.maxBlockSize = 1000;

}

#endif
