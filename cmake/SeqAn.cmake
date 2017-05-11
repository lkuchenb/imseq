# Check if the user specified a location
if (SEQAN_ROOT)
	if (NOT EXISTS ${SEQAN_ROOT}/util/cmake/seqan-config.cmake
			AND NOT EXISTS ${SEQAN_ROOT}/lib/cmake/seqan/seqan-config.cmake )
		message ( FATAL_ERROR "Invalid SEQAN_ROOT '${SEQAN_ROOT}'." )
	endif()
else()
	set ( SEQAN_URL "https://github.com/seqan/seqan/releases/download/seqan-v2.3.2/seqan-library-2.3.2.zip")
	set ( SEQAN_MD5 "d92e2d813b4939332a9f64c97be399b9")
	set ( SEQAN_ZIP_OUT ${CMAKE_CURRENT_BINARY_DIR}/seqan-library-2.3.2.zip )
	set ( SEQAN_ROOT ${CMAKE_CURRENT_BINARY_DIR}/seqan-library-2.3.2 )
	
	message ( STATUS "SEQAN_ROOT set to '${SEQAN_ROOT}'" )

	# Check if the desired version was already downloaded and unpacked. If
	# not, perform the download and unpacking.
	if (NOT EXISTS ${SEQAN_ROOT}/lib/cmake/seqan/seqan-config.cmake )
		# Download zip file
		message ( STATUS "Downloading ${SEQAN_URL}" )
		file (DOWNLOAD ${SEQAN_URL} ${SEQAN_ZIP_OUT}
			EXPECTED_MD5 ${SEQAN_MD5}
			SHOW_PROGRESS STATUS status)
		list ( GET status 0 ret )
		list ( GET status 0 str)
		if ( NOT ret EQUAL 0)
			message (FATAL_ERROR "Download failed" )
		endif()
		# Unpack zip file
		message ( STATUS "Unpacking ${SEQAN_ZIP_OUT}")
		execute_process(COMMAND cmake -E tar zxf ${SEQAN_ZIP_OUT}  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
		# Remove zip file
		if ( EXISTS ${SEQAN_ZIP_OUT} )
			file (REMOVE ${SEQAN_ZIP_OUT})
		endif()
	endif()
	# Check if FindSeqAn.cmake can be found wher it should
	if (NOT EXISTS ${SEQAN_ROOT}/lib/cmake/seqan/seqan-config.cmake )
		message (FATAL_ERROR "Failed to download and unpack '${SEQAN_URL}'")
	endif()
endif ()

# Check whether SEQAN_ROOT root points to the repository version or the
# library-only version of seqan
if ( EXISTS ${SEQAN_ROOT}/util/cmake/seqan-config.cmake )
	# seqan_root points to a repository clone or copy
	set ( SeqAn_DIR ${SEQAN_ROOT}/util/cmake/ )
elseif ( EXISTS ${SEQAN_ROOT}/lib/cmake/seqan/seqan-config.cmake )
	# seqan_root points to a library-only folder
	set ( SeqAn_DIR ${SEQAN_ROOT}/lib/cmake/seqan/ )
endif()

message ( STATUS "SeqAn_DIR='${SeqAn_DIR}'" )

# Make sure the SeqAn module configures the correct SeqAn installation
SET( SEQAN_INCLUDE_PATH ${SEQAN_ROOT}/include/ )

# Load the SeqAn module and fail if not found.
set ( SEQAN_DISABLE_VERSION_CHECK )
find_package (SeqAn REQUIRED)
