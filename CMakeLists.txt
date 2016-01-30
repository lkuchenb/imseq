# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/imseq/CMakeLists.txt
#
# CMakeLists.txt file for imseq.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_apps_imseq)
message (STATUS "Configuring apps/imseq")

# Enable C++11 support if possible
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    else ()
        message(FATAL_ERROR "To compile with GCC, GCC 4.7 or newer is required! You can set the C++ compiler using '-DCMAKE_CXX_COMPILER=/path/to/compiler' when invoking cmake!")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
else()
    message(FATAL_ERROR "Only the following compilers are supported and have to be C++11 compliant: gcc, clang. You can set the C++ compiler using '-DCMAKE_CXX_COMPILER=/path/to/compiler' when invoking cmake!")
endif ()

# ----------------------------------------------------------------------------
# Add subdirectories
# ----------------------------------------------------------------------------

add_subdirectory(src)

add_subdirectory(unit_test)
