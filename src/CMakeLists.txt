cmake_minimum_required (VERSION 3.0.0)
project (imseq)
message ( STATUS "Configuring IMSEQ" )


# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

find_package ( ZLIB REQUIRED )
find_package (Git)
include ( SeqAn )

if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --dirty=-dirty
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE DESCRIBE_BUILD
	ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    #    string(REGEX MATCH "[0-9]+" BUILD_NUMBER ${DESCRIBE_BUILD})
endif()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Disable SEQAN version check unless explicitly requested otherwise
if ( NOT ALLOW_SEQAN_VERSION_CHECK )
	add_definitions ( -DSEQAN_DISABLE_VERSION_CHECK )
endif()

# IMSEQ static build for Linux
if ("${CMAKE_IMSEQ_LINUX_STATIC}" EQUAL "true")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    SET(BUILD_SHARED_LIBRARIES OFF)
    SET(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
endif()

# IMSEQ executable
add_executable (imseq 
	aa_translate.h
	barcode_correction.h
	cdr3_cli.h
	cdr_utils.h
	clone.h
	clone_store.h
	cluster_log.cpp
	cluster_log.h
	cluster_result.cpp
	cluster_result.h
	collection_utils.h
	extdir_oldir_conversion.h
	fastq_io.h
	fastq_io_types.h
	fastq_multi_record.h
	fastq_multi_record_types.h
	file_utils.h
	fixed_size_types.h
	globalData.h
	imseq.cpp
	imseq.h
	logging.cpp
	logging.h
	overlap_specs.h
	progress_bar.cpp
	progress_bar.h
	qc_basics.h
	referencePreparation.h
	reject.h
	runtime_options.h
	segment_ambiguity.h
	segment_meta.cpp
	segment_meta.h
	sequence_data.h
	sequence_data_types.h
	thread_check.h
	thread_pool.cpp
	thread_pool.h
	timeFormat.h
	version_number.cpp
	version_number.h
	vjMatching.h
	)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (imseq ${SEQAN_LIBRARIES})

install (TARGETS imseq
         DESTINATION bin)

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "-Wno-address ${CMAKE_CXX_FLAGS}")
