# FindTargetEigen3.cmake
# --------------------
#
# Eigen3 cmake module to wrap Eigen3 suitable for Libint2, whether Eigen3Config, FindEigen3, or raw variables, in a target.
#
# This module sets the following variables in your project: ::
#
##   TargetEigen3_FOUND - true if Eigen3 and all required components found on the system
##   TargetEigen3_VERSION - Eigen3 version in format Major.Minor.Release
##   TargetEigen3_MESSAGE - status message with Eigen3 library path list and version
##
## Note that components are passed along to find_package(HDF5 (untested) but not checked in the direct TargetHDF5Config
## Note that version checking/attaching not working yet
##
## This module *unsets* the following conventional HDF5 variables so as
##   to force using the target: ::
##
##   HDF5_FOUND
##   HDF5_VERSION
##   HDF5_INCLUDE_DIRS
##   HDF5_LIBRARIES
##
## Exported targets::
##
## If Eigen3 is found, this module defines the following :prop_tgt:`IMPORTED`
## target. ::
##
##   tgt::Eigen3 - the Eigen3 libraries with headers attached.
##
## Suggested usage::
##
##   find_package(TargetHDF5)
##   find_package(TargetHDF5 1.8.16 REQUIRED)
#
#
# The following variables can be set to guide the search for this package::
#
#   TargetEigen3_DIR - CMake variable, set to directory containing this Config file
#   CMAKE_PREFIX_PATH - CMake variable, set to root directory of this package

set(PN TargetEigen3)

# 1st precedence - libraries passed in through -DHDF5_LIBRARIES
if (HDF5_LIBRARIES AND HDF5_INCLUDE_DIRS)
    if (HDF5_VERSION)
        if (NOT ${PN}_FIND_QUIETLY)
            message (STATUS "HDF5 detection suppressed.")
        endif()

        add_library (tgt::hdf5 INTERFACE IMPORTED)
        set_property (TARGET tgt::hdf5 PROPERTY INTERFACE_LINK_LIBRARIES ${HDF5_LIBRARIES})
        set_property (TARGET tgt::hdf5 PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
        set (${PN}_VERSION ${HDF5_VERSION})
    else()
        message (FATAL_ERROR "Humor the build system - pass in the version, too (for example, -DHDF5_VERSION=1.8.17).")
    endif()
else()
    # 2nd precedence - target already prepared and findable in TargetEigen3Config.cmake
    find_package (TargetEigen3 QUIET CONFIG)
    if ((TARGET tgt::Eigen) AND (${PN}_VERSION))
        if (NOT ${PN}_FIND_QUIETLY)
            message (STATUS "TargetEigen3Config detected.")
        endif()
    else()
        # 3rd precedence - usual variables from Eigen3Config.cmake
        find_package (Eigen3 QUIET CONFIG)
#        find_package (Eigen3 QUIET MODULE)
# Eigen3::Eigen
        if (NOT ${PN}_FIND_QUIETLY)
            message (STATUS "Eigen3 detected.")
        endif()
    
        add_library (tgt::Eigen INTERFACE IMPORTED)

        get_property(_iid1 TARGET Eigen3::Eigen PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
        set_property (TARGET tgt::Eigen PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${_iid1})
# set_property (TARGET tgt::eigen PROPERTY INTERFACE_LINK_LIBRARIES ${HDF5_LIBRARIES})
        #set_property (TARGET tgt::hdf5 PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
        set (${PN}_VERSION ${Eigen3_VERSION})

        #unset (HDF5_FOUND)
        #unset (HDF5_VERSION)
        #unset (HDF5_LIBRARIES)
        #unset (HDF5_INCLUDE_DIRS)
    endif()
endif()    

get_property(_iid TARGET tgt::Eigen PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
set(${PN}_MESSAGE "Found Eigen3: ${_iid} (found version ${${PN}_VERSION})")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(${PN}
                                  REQUIRED_VARS ${PN}_MESSAGE
                                  VERSION_VAR ${PN}_VERSION)

# FindTargetLAPACK.cmake
# ----------------------
#
# LAPACK cmake module to wrap FindLAPACK.cmake in a target.
#
# This module sets the following variables in your project: ::
#
#   TargetLAPACK_FOUND - true if BLAS/LAPACK found on the system
#   TargetLAPACK_MESSAGE - status message with BLAS/LAPACK library path list
#
# This module *unsets* the following conventional LAPACK variables so as
#   to force using the target: ::
#
#   LAPACK_FOUND
#   LAPACK_LIBRARIES
#
# In order of decreasing precedence, this module returns in a target ``tgt::lapack``
#  (1) the libraries passed through CMake variable LAPACK_LIBRARIES,
#  (2) the libraries defined in a detectable TargetLAPACKConfig.cmake file
#      (skip via DISABLE_FIND_PACKAGE_TargetLAPACK), or
#  (3) the libraries detected by the usual FindLAPACK.cmake module.
#

#set(PN TargetLAPACK)
#
## 1st precedence - libraries passed in through -DLAPACK_LIBRARIES
#if (LAPACK_LIBRARIES)
#    if (NOT ${PN}_FIND_QUIETLY)
#        message (STATUS "LAPACK detection suppressed.")
#    endif()
#
#    add_library (tgt::lapack INTERFACE IMPORTED)
#    set_property (TARGET tgt::lapack PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
#else()
#    # 2nd precedence - target already prepared and findable in TargetLAPACKConfig.cmake
#    if (NOT "${DISABLE_FIND_PACKAGE_${PN}}")
#        find_package (TargetLAPACK QUIET CONFIG)
#    endif()
#    if (TARGET tgt::lapack)
#        if (NOT ${PN}_FIND_QUIETLY)
#            message (STATUS "TargetLAPACKConfig detected.")
#        endif()
#    else()
#        # 3rd precedence - usual variables from FindLAPACK.cmake
#        find_package (LAPACK QUIET MODULE)
#        if (NOT ${PN}_FIND_QUIETLY)
#            message (STATUS "LAPACK detected.")
#        endif()
#    
#        add_library (tgt::lapack INTERFACE IMPORTED)
#        set_property (TARGET tgt::lapack PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
#
#        unset (LAPACK_FOUND)
#        unset (LAPACK_LIBRARIES)
#    endif()
#endif()    
#
#get_property (_ill TARGET tgt::lapack PROPERTY INTERFACE_LINK_LIBRARIES)
#set (${PN}_MESSAGE "Found LAPACK: ${_ill}")
#if ((TARGET tgt::blas) AND (TARGET tgt::lapk))
#    get_property (_illb TARGET tgt::blas PROPERTY INTERFACE_LINK_LIBRARIES)
#    get_property (_illl TARGET tgt::lapk PROPERTY INTERFACE_LINK_LIBRARIES)
#    set (${PN}_MESSAGE "Found LAPACK: ${_illl};${_illb}")
#endif()
#
#include(FindPackageHandleStandardArgs)
#find_package_handle_standard_args (${PN} DEFAULT_MSG ${PN}_MESSAGE)
#
