cmake_minimum_required (VERSION 3.0.0)

set ( SEQAN_MOD_DIR "${SEQAN_ROOT}/util/cmake/" )
if ( EXISTS "${SEQAN_MOD_DIR}" )
	message ( STATUS "Configuring IMSEQ unit tests" )
	LIST ( APPEND CMAKE_MODULE_PATH "${SEQAN_MOD_DIR}" )
	include (SeqAnBuildSystem)

	seqan_setup_library ()
	seqan_register_tests ()

	include_directories (${SEQAN_INCLUDE_DIRS})
	add_definitions (${SEQAN_DEFINITIONS})
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

	# unit_tests_imseq executable
	add_executable (unit_tests_imseq
		unit_tests_imseq.cpp
		unit_tests_imseq_barcode_correction.h
		unit_tests_imseq_fastq_io.h
		unit_tests_imseq_fastq_multi_record.h
		unit_tests_imseq_qc_basics.h
		)

	# Add dependencies found by find_package (SeqAn).
	target_link_libraries (unit_tests_imseq ${SEQAN_LIBRARIES})
else()
	message ( STATUS "Not configuring unit tests. Set SEQAN_ROOT to a working copy of the SeqAn repository if needed." )
endif()


