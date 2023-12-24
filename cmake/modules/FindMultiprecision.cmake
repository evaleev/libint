# FindMultiprecision.cmake
# ------------------------
#
# GMP, GMP++, MPIR, MPFR CMake module for Unix and Windows.
#
# Built upon FindMPFR.cmake by
#   Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
#   Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
#   Copyright (c) 2010 Jitse Niesen, <jitse@maths.leeds.ac.uk>
#   Copyright (c) 2015 Jack Poulson, <jack.poulson@gmail.com>
#   Redistribution and use is allowed according to the terms of the BSD license.
#
#
# This module sets the following variables in your project:
#
# ::
#
#   Multiprecision_FOUND - true if all required components of MPFR found on the system with correct version
#   MPFR_VERSION - MPFR version if mpfr library detected
#
#
# Available components:
#
# ::
#
#   gmp - search for at least Multiprecision::gmp target
#   gmpxx - search for at least Multiprecision::gmpxx target
#   mpfr - search for at least Multiprecision::mpfr target
#
#
# Exported targets:
#
# ::
#
# If no components are requested, this module defines at least the following
# :prop_tgt:`IMPORTED` target. ::
#
#   Multiprecision::gmp - gmp.h and GMP library
#
# Depending on components requested, this module defines up to the following
# :prop_tgt:`IMPORTED` targets. ::
#
#   Multiprecision::gmp - gmp.h and C GMP library
#   Multiprecision::mpfr - mpfr.h and MPFR library (needs Multiprecision::gmp)
#   Multiprecision::gmpxx - gmpxx.h and C++ GMP library (needs Multiprecision::gmp)
#
#
# Suggested usage:
#
# ::
#
#   find_package(Multiprecision)
#   find_package(Multiprecision 4.0.0 REQUIRED COMPONENTS mpfr)
#
#
# The following variables can be set to guide the search for this package:
#
# ::
#
#   Multiprecision_ROOT - CMake variable, set to root directory of this package
#   CMAKE_PREFIX_PATH - CMake variable, set to root directory of this package

# if no components given, act like a FindGMP
list(LENGTH Multiprecision_FIND_COMPONENTS _lcomp)
if (_lcomp EQUAL 0)
    list(APPEND Multiprecision_FIND_COMPONENTS gmp)
endif()

if(WIN32)
    list(INSERT CMAKE_FIND_LIBRARY_SUFFIXES 0 ".dll")
endif()

find_path(
  MPFR_INCLUDE
  NAMES
    mpfr.h
  PATHS
    $ENV{GMPDIR}
    $ENV{MPFRDIR}
    ${INCLUDE_INSTALL_DIR}
  )

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

find_library(
  MPFR_LIBRARY
  NAMES
    mpfr
  PATHS
    $ENV{GMPDIR}
    $ENV{MPFRDIR}
    ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
    bin
    ${MPFR_ROOT}/bin
  )

find_path(
  GMP_INCLUDE
  NAMES
    gmp.h
  PATHS
    $ENV{GMPDIR}
    $ENV{MPFRDIR}
    ${INCLUDE_INSTALL_DIR}
  )

find_library(
  GMP_LIBRARY
  NAMES
    gmp
  PATHS
    $ENV{GMPDIR}
    $ENV{MPFRDIR}
    ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
    bin
    ${MPFR_ROOT}/bin
  )

find_path(
  GMPXX_INCLUDE
  NAMES
    gmpxx.h
  PATHS
    $ENV{GMPDIR}
    $ENV{MPFRDIR}
    ${INCLUDE_INSTALL_DIR}
  )

find_library(
  GMPXX_LIBRARY
  NAMES
    gmpxx
    gmp  # gmp.dll on Win c-f conda package contains cxx (actually a copy of mpir, a drop-in replacement for gmp)
  PATHS
    $ENV{GMPDIR}
    $ENV{MPFRDIR}
    ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
    bin
    ${MPFR_ROOT}/bin
  )

if (CMAKE_VERSION VERSION_GREATER 3.15)
    message(VERBOSE "GMP    inc=${GMP_INCLUDE}   lib=${GMP_LIBRARY}")
    message(VERBOSE "GMPXX  inc=${GMPXX_INCLUDE} lib=${GMPXX_LIBRARY}")
    message(VERBOSE "MPFR   inc=${MPFR_INCLUDE}  lib=${MPFR_LIBRARY}")
endif()

if (GMP_INCLUDE AND GMP_LIBRARY)
    set(Multiprecision_gmp_FOUND 1)  # Multiprecision::gmp

    if (MPFR_INCLUDE AND MPFR_LIBRARY AND MPFR_VERSION_OK)
        set(Multiprecision_mpfr_FOUND 1)  # Multiprecision::mpfr
    endif()

    if (GMPXX_INCLUDE AND GMPXX_LIBRARY)
        set(Multiprecision_gmpxx_FOUND 1)  # Multiprecision::gmpxx
    endif()
endif()

# thanks, https://stackoverflow.com/a/9328525
function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    set(founds "")
    foreach (_variableName ${_variableNames})
        if (ARGV0)
            unset(MATCHED)
            string(REGEX MATCH ${ARGV0} MATCHED ${_variableName})
            if (NOT MATCHED)
                continue()
            endif()
            if (NOT ${${_variableName}})
                continue()
            endif()
        endif()
        list(APPEND founds ${CMAKE_MATCH_1})
    endforeach()
    message(STATUS "${ARGV1}${founds}")
endfunction()

macro(check_required_components _NAME)
    foreach(comp ${${_NAME}_FIND_COMPONENTS})
        if(NOT ${_NAME}_${comp}_FOUND)
            if(${_NAME}_FIND_REQUIRED_${comp})
                set(${_NAME}_FOUND FALSE)
            endif()
        endif()
    endforeach()
endmacro()

if(NOT CMAKE_REQUIRED_QUIET)
    message(STATUS "FindMultiprecision components requested: ${Multiprecision_FIND_COMPONENTS}")
    dump_cmake_variables("^Multiprecision_([A-Za-z0-9_]+)_FOUND$" "FindMultiprecision components found: ")
endif()

include(FindPackageHandleStandardArgs)
set(Multiprecision_FOUND 1)
check_required_components(Multiprecision)
find_package_handle_standard_args(
  Multiprecision
  REQUIRED_VARS  # can be removed upon CMake 3.18
    GMP_LIBRARY
    GMP_INCLUDE
    Multiprecision_FOUND
  VERSION_VAR MPFR_VERSION
  HANDLE_COMPONENTS
  )

if(WIN32)
  string(REPLACE ".lib" ".dll" GMP_LIBRARY_DLL "${GMP_LIBRARY}")
  string(REPLACE ".lib" ".dll" GMPXX_LIBRARY_DLL "${GMPXX_LIBRARY}")
  string(REPLACE ".lib" ".dll" MPFR_LIBRARY_DLL "${MPFR_LIBRARY}")
endif()

# now that `find_package(Multiprecision COMPONENTS ...)` will succeed, create targets
if ((Multiprecision_gmp_FOUND EQUAL 1) AND NOT TARGET Multiprecision::gmp)
    if (EXISTS "${GMP_LIBRARY_DLL}")
        add_library(Multiprecision::gmp SHARED IMPORTED)
        set_target_properties(
          Multiprecision::gmp
          PROPERTIES
            IMPORTED_LOCATION "${GMP_LIBRARY_DLL}"
            IMPORTED_IMPLIB "${GMP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE}"
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
          )
    else()
        add_library(Multiprecision::gmp UNKNOWN IMPORTED)
        set_target_properties(
          Multiprecision::gmp
          PROPERTIES
            IMPORTED_LOCATION "${GMP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE}"
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
          )
    endif()
endif()

if ((Multiprecision_mpfr_FOUND EQUAL 1) AND NOT TARGET Multiprecision::mpfr)
    if (EXISTS "${MPFR_LIBRARY_DLL}")
        add_library(Multiprecision::mpfr SHARED IMPORTED)
        set_target_properties(
          Multiprecision::mpfr
          PROPERTIES
            IMPORTED_LOCATION "${MPFR_LIBRARY_DLL}"
            IMPORTED_IMPLIB "${MPFR_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE}"
            INTERFACE_LINK_LIBRARIES Multiprecision::gmp
          )
    else()
        add_library(Multiprecision::mpfr UNKNOWN IMPORTED)
        set_target_properties(
          Multiprecision::mpfr
          PROPERTIES
            IMPORTED_LOCATION "${MPFR_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE}"
            INTERFACE_LINK_LIBRARIES Multiprecision::gmp
          )
    endif()
endif()

if ((Multiprecision_gmpxx_FOUND EQUAL 1) AND NOT TARGET Multiprecision::gmpxx)
    if (EXISTS "${GMPXX_LIBRARY_DLL}")
        add_library(Multiprecision::gmpxx SHARED IMPORTED)
        set_target_properties(
          Multiprecision::gmpxx
          PROPERTIES
            IMPORTED_LOCATION "${GMPXX_LIBRARY_DLL}"
            IMPORTED_IMPLIB "${GMPXX_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMPXX_INCLUDE}"
            IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
            INTERFACE_LINK_LIBRARIES Multiprecision::gmp
          )
    else()
        add_library(Multiprecision::gmpxx UNKNOWN IMPORTED)
        set_target_properties(
          Multiprecision::gmpxx
          PROPERTIES
            IMPORTED_LOCATION "${GMPXX_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMPXX_INCLUDE}"
            IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
            INTERFACE_LINK_LIBRARIES Multiprecision::gmp
          )
    endif()
endif()

