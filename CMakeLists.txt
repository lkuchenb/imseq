cmake_minimum_required (VERSION 3.0.0)

# Add custom module folder
LIST ( APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/ )

# Find or download SeqAn
include ( SeqAn )

# Make sure executables end up directly in the build folder
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# IMSEQ
add_subdirectory(src)

# Unit tests
add_subdirectory(unit_test)
