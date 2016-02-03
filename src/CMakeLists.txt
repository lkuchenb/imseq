cmake_minimum_required (VERSION 2.8.2)
project (seqan_apps_imseq_src)
message (STATUS "Configuring apps/imseq/src")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)
find_package (Git)

if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --dirty=-dirty
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE DESCRIBE_BUILD
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    #    string(REGEX MATCH "[0-9]+" BUILD_NUMBER ${DESCRIBE_BUILD})
endif()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})
add_definitions (-DBUILD_NUMBER=${DESCRIBE_BUILD})

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
	fastq_multi_record_types.cpp
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

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (imseq)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install imseq in ${PREFIX}/bin directory
install (TARGETS imseq
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/imseq for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (imseq)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:imseq")
  set (CPACK_PACKAGE_NAME "imseq")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "imseq")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (imseq "imseq")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:imseq")

