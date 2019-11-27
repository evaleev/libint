# Try to find the MPFR library
# See http://www.mpfr.org/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(MPFR 2.3.0)
# to require version 2.3.0 to newer of MPFR.
#
# Once done this will define
#
#  MPFR_FOUND - system has MPFR lib with correct version
#  MPFR_INCLUDE - the MPFR include directory
#  MPFR_LIBRARY - the MPFR library
#  MPFR_VERSION - MPFR version
#  MPFR::GMP
#  MPFR::GMPXX
#  MPFR::MPFR
#  MPFR::MPFRXX

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2010 Jitse Niesen, <jitse@maths.leeds.ac.uk>
# Copyright (c) 2015 Jack Poulson, <jack.poulson@gmail.com>
# Redistribution and use is allowed according to the terms of the BSD license.

message("prefix: ${CMAKE_FIND_LIBRARY_PREFIXES}")
message("suffix: ${CMAKE_FIND_LIBRARY_SUFFIXES}")
list(APPEND CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
message("suffix: ${CMAKE_FIND_LIBRARY_SUFFIXES}")
message("root: ${MPFR_ROOT}")

find_path(MPFR_INCLUDE NAMES mpfr.h PATHS $ENV{GMPDIR} $ENV{MPFRDIR}
        ${INCLUDE_INSTALL_DIR})

# Set MPFR_FIND_VERSION to 1.0.0 if no minimum version is specified
if(NOT MPFR_FIND_VERSION)
    if(NOT MPFR_FIND_VERSION_MAJOR)
        set(MPFR_FIND_VERSION_MAJOR 1)
    endif()
    if(NOT MPFR_FIND_VERSION_MINOR)
        set(MPFR_FIND_VERSION_MINOR 0)
    endif()
    if(NOT MPFR_FIND_VERSION_PATCH)
        set(MPFR_FIND_VERSION_PATCH 0)
    endif()
    set(MPFR_FIND_VERSION
            "${MPFR_FIND_VERSION_MAJOR}.${MPFR_FIND_VERSION_MINOR}.${MPFR_FIND_VERSION_PATCH}")
endif()

if(MPFR_INCLUDE)
    # Query MPFR_VERSION
    file(READ "${MPFR_INCLUDE}/mpfr.h" _mpfr_version_header)

    string(REGEX MATCH "define[ \t]+MPFR_VERSION_MAJOR[ \t]+([0-9]+)"
            _mpfr_major_version_match "${_mpfr_version_header}")
    set(MPFR_MAJOR_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+MPFR_VERSION_MINOR[ \t]+([0-9]+)"
            _mpfr_minor_version_match "${_mpfr_version_header}")
    set(MPFR_MINOR_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+MPFR_VERSION_PATCHLEVEL[ \t]+([0-9]+)"
            _mpfr_patchlevel_version_match "${_mpfr_version_header}")
    set(MPFR_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")

    set(MPFR_VERSION
            ${MPFR_MAJOR_VERSION}.${MPFR_MINOR_VERSION}.${MPFR_PATCHLEVEL_VERSION})

    # Check whether found version exceeds minimum required
    if(${MPFR_VERSION} VERSION_LESS ${MPFR_FIND_VERSION})
        set(MPFR_VERSION_OK FALSE)
        message(STATUS "MPFR version ${MPFR_VERSION} found in ${MPFR_INCLUDE}, "
                "but at least version ${MPFR_FIND_VERSION} is required")
    else()
        set(MPFR_VERSION_OK TRUE)
    endif()
endif()

find_library(MPFR_LIBRARY
    NAMES mpfr
#          mpfr_static
    PATHS $ENV{GMPDIR} $ENV{MPFRDIR} ${LIB_INSTALL_DIR}
    PATH_SUFFIXES bin ${MPFR_ROOT}/bin)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG
        MPFR_INCLUDE MPFR_LIBRARY MPFR_VERSION_OK)
mark_as_advanced(MPFR_INCLUDE MPFR_LIBRARY)

if(MPFR_INCLUDE AND MPFR_LIBRARY AND NOT TARGET MPFR::Library)
    add_library(MPFR::Library INTERFACE IMPORTED)
    set_target_properties(MPFR::Library
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${MPFR_INCLUDE}
                   INTERFACE_LINK_LIBRARIES ${MPFR_LIBRARY})
endif()

# from GMP, Libint2 build_libint needs the C++ header and MPFR needs the C library
#                           test needs the C++ library
find_path(GMP_INCLUDE
    NAMES gmp.h
    PATHS $ENV{GMPDIR} $ENV{MPFRDIR} ${INCLUDE_INSTALL_DIR})
find_library(GMP_LIBRARY
    NAMES gmp
#          gmp_static
    PATHS $ENV{GMPDIR} $ENV{MPFRDIR} ${LIB_INSTALL_DIR}
    PATH_SUFFIXES bin ${MPFR_ROOT}/bin)

if (GMP_INCLUDE AND GMP_LIBRARY AND NOT TARGET MPFR::GMP)
    add_library(MPFR::GMP INTERFACE IMPORTED)
    set_target_properties(MPFR::GMP
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDE}
                   INTERFACE_LINK_LIBRARIES ${GMP_LIBRARY})
endif()

if (TARGET MPFR::GMP AND TARGET MPFR::Library AND NOT TARGET MPFR::MPFR)
    add_library(MPFR::MPFR INTERFACE IMPORTED)
    target_link_libraries(MPFR::MPFR INTERFACE MPFR::Library MPFR::GMP)
endif()


find_path(GMPXX_INCLUDE
    NAMES gmpxx.h
    PATHS $ENV{GMPDIR} $ENV{MPFRDIR} ${INCLUDE_INSTALL_DIR})
find_library(GMPXX_LIBRARY
    NAMES gmpxx
#          gmpxx_static
    PATHS $ENV{GMPDIR} $ENV{MPFRDIR} ${LIB_INSTALL_DIR}
    PATH_SUFFIXES bin ${MPFR_ROOT}/bin)

if (GMPXX_INCLUDE AND GMPXX_LIBRARY AND TARGET MPFR::GMP AND NOT TARGET MPFR::GMPXX)
    add_library(MPFR::GMPXX INTERFACE IMPORTED)
    set_target_properties(MPFR::GMPXX
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${GMPXX_INCLUDE}
                   INTERFACE_LINK_LIBRARIES ${GMPXX_LIBRARY})
    target_link_libraries(MPFR::GMPXX INTERFACE MPFR::GMP)
endif()

if (TARGET MPFR::GMPXX AND TARGET MPFR::Library AND NOT TARGET MPFR::MPFRXX)
    add_library(MPFR::MPFRXX INTERFACE IMPORTED)
    target_link_libraries(MPFR::MPFRXX INTERFACE MPFR::Library MPFR::GMPXX)
endif()
