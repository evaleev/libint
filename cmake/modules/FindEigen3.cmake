# FindEigen3.cmake
# ----------------
#
# Eigen3 cmake module to wrap Eigen3 suitable for Libint2. Copied from Eigen v3.4.0 source and modified as follows:
# * Added `NO_CMAKE_PACKAGE_REGISTRY` to `find_package(Eigen3 ...)` to avoid issues with wiped build
#   directory when looking for installed eigen. Eigen3 registers its *build* tree with the user package registry.
# * Added `LIBINT2_LOCAL_Eigen3_FIND` block to forcibly load hard-coded Eigen3 location detected during Libint2 library build.
# * Move default Eigen3_FIND_VERSION_* from 2.91.0 to 3.0.0 so that it doesn't reject Eigen v3 installations.
#
# - Try to find Eigen3 lib
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Eigen3 3.1.2)
# to require version 3.1.2 or newer of Eigen3.
#
# Once done this will define
#
#  EIGEN3_FOUND - system has eigen lib with correct version
#  EIGEN3_INCLUDE_DIR - the eigen include directory
#  EIGEN3_VERSION - eigen version
#
# and the following imported target:
#
#  Eigen3::Eigen - The header-only Eigen library
#
# This module reads hints about search locations from 
# the following environment variables:
#
# EIGEN3_ROOT
# EIGEN3_ROOT_DIR

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.

if(NOT Eigen3_FIND_VERSION)
  if(NOT Eigen3_FIND_VERSION_MAJOR)
    set(Eigen3_FIND_VERSION_MAJOR 3)
  endif()
  if(NOT Eigen3_FIND_VERSION_MINOR)
    set(Eigen3_FIND_VERSION_MINOR 0)
  endif()
  if(NOT Eigen3_FIND_VERSION_PATCH)
    set(Eigen3_FIND_VERSION_PATCH 0)
  endif()

  set(Eigen3_FIND_VERSION "${Eigen3_FIND_VERSION_MAJOR}.${Eigen3_FIND_VERSION_MINOR}.${Eigen3_FIND_VERSION_PATCH}")
endif()

macro(_eigen3_check_version)
  file(READ "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" _eigen3_version_header)

  string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _eigen3_world_version_match "${_eigen3_version_header}")
  set(EIGEN3_WORLD_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen3_major_version_match "${_eigen3_version_header}")
  set(EIGEN3_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen3_minor_version_match "${_eigen3_version_header}")
  set(EIGEN3_MINOR_VERSION "${CMAKE_MATCH_1}")

  set(EIGEN3_VERSION ${EIGEN3_WORLD_VERSION}.${EIGEN3_MAJOR_VERSION}.${EIGEN3_MINOR_VERSION})
  if(${EIGEN3_VERSION} VERSION_LESS ${Eigen3_FIND_VERSION})
    set(EIGEN3_VERSION_OK FALSE)
  else()
    set(EIGEN3_VERSION_OK TRUE)
  endif()

  if(NOT EIGEN3_VERSION_OK)

    message(STATUS "Eigen3 version ${EIGEN3_VERSION} found in ${EIGEN3_INCLUDE_DIR}, "
                   "but at least version ${Eigen3_FIND_VERSION} is required")
  endif()
endmacro()

if (LIBINT2_LOCAL_Eigen3_FIND)
    include("${CMAKE_CURRENT_LIST_DIR}/libint2-targets-eigen3.cmake")

    if (TARGET Libint2::Eigen)
        get_property(_EIGEN3_INCLUDE_DIRS TARGET Libint2::Eigen PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
        list(GET _EIGEN3_INCLUDE_DIRS 0 EIGEN3_INCLUDE_DIR)

        _eigen3_check_version()
        set(EIGEN3_FOUND ${EIGEN3_VERSION_OK})
        set(Eigen3_FOUND ${EIGEN3_VERSION_OK})

        add_library(Eigen3::Eigen ALIAS Libint2::Eigen)
    else()
        message(STATUS "Eigen3 exact installation detected/used by Libint library build requested "
                       "from ${CMAKE_CURRENT_LIST_DIR}/libint2-targets-eigen3.cmake but failed.")
    endif()

elseif (EIGEN3_INCLUDE_DIR)

  # in cache already
  _eigen3_check_version()
  set(EIGEN3_FOUND ${EIGEN3_VERSION_OK})
  set(Eigen3_FOUND ${EIGEN3_VERSION_OK})

else ()
  
  # search first if an Eigen3Config.cmake is available in the system,
  # if successful this would set EIGEN3_INCLUDE_DIR and the rest of
  # the script will work as usual
  find_package(Eigen3 ${Eigen3_FIND_VERSION} NO_MODULE QUIET NO_CMAKE_PACKAGE_REGISTRY)

  if(NOT EIGEN3_INCLUDE_DIR)
    find_path(EIGEN3_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
        HINTS
        ENV EIGEN3_ROOT 
        ENV EIGEN3_ROOT_DIR
        PATHS
        ${CMAKE_INSTALL_PREFIX}/include
        ${KDE4_INCLUDE_DIR}
        PATH_SUFFIXES eigen3 eigen
      )
  endif()

  if(EIGEN3_INCLUDE_DIR)
    _eigen3_check_version()
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Eigen3 DEFAULT_MSG EIGEN3_INCLUDE_DIR EIGEN3_VERSION_OK)

  mark_as_advanced(EIGEN3_INCLUDE_DIR)

endif()

if(EIGEN3_FOUND AND NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
endif()

