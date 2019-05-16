# FindLibint2.cmake
#
# Finds the Libint2 header and library using the pkg-config program.
#
# Use this module by invoking find_package() as follows:
#   find_package(Libint2
#                [version] [EXACT]      # Minimum or EXACT version e.g. 2.5.0
#                [REQUIRED]             # Fail with error if Libint2 is not found
#               )
#
# The behavior can be controlled by setting the following variables
#
#    LIBINT2_SHARED_LIBRARY_ONLY              if true, will look for shared lib only; may be needed for some platforms
#                                             where linking errors or worse, e.g. duplicate static data, occur if
#                                             linking shared libraries against static Libint2 library.
#    PKG_CONFIG_PATH (environment variable)   Add the libint2 install prefix directory (e.g. /usr/local)
#                                             to specify where to look for libint2
#    CMAKE_MODULE_PATH                        Add the libint2 install prefix directory (e.g. /usr/local)
#                                             to specify where to look for libint2
#
# This will define the following CMake cache variables
#
#    LIBINT2_FOUND           - true if libint2.h header and libint2 library were found
#    LIBINT2_VERSION         - the libint2 version
#    LIBINT2_INCLUDE_DIRS    - (deprecated: use the CMake IMPORTED targets listed below) list of libint2 include directories
#    LIBINT2_LIBRARIES       - (deprecated: use the CMake IMPORTED targets listed below) list of libint2 libraries
#
# and the following imported targets
#
#     Libint2::int2          - library only
#     Libint2::cxx           - library + C++11 API; may need to add dependence on Eigen3 and/or Boost.Preprocessor if
#                              was not found by Libint at configure time
#
# Author: Eduard Valeyev - libint@valeyev.net

# need cmake 3.8 for cxx_std_11 compile feature
if(CMAKE_VERSION VERSION_LESS 3.8.0)
    message(FATAL_ERROR "This file relies on consumers using CMake 3.8.0 or greater.")
endif()

find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
    # point pkg-config to the location of this tree's libint2.pc
    set(ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_LIST_DIR}/../../pkgconfig:$ENV{PKG_CONFIG_PATH}")
    if(LIBINT2_FIND_QUIETLY)
        pkg_check_modules(PC_LIBINT2 QUIET libint2)
    else()
        pkg_check_modules(PC_LIBINT2 libint2)
    endif()
    set(LIBINT2_VERSION ${PC_LIBINT2_VERSION})

    find_path(LIBINT2_INCLUDE_DIR
            NAMES libint2.h
            PATHS ${PC_LIBINT2_INCLUDE_DIRS}
            PATH_SUFFIXES libint2
            )

    if (LIBINT2_SHARED_LIBRARY_ONLY)
        set(_LIBINT2_LIB_NAMES "libint2.so" "libint2.dylib")
    else (LIBINT2_SHARED_LIBRARY_ONLY)
        set(_LIBINT2_LIB_NAMES "int2")
    endif(LIBINT2_SHARED_LIBRARY_ONLY)

    find_library(LIBINT2_LIBRARY NAMES ${_LIBINT2_LIB_NAMES} HINTS ${PC_LIBINT2_LIBRARY_DIRS})

    mark_as_advanced(LIBINT2_FOUND LIBINT2_INCLUDE_DIR LIBINT2_LIBRARY LIBINT2_VERSION)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Libint2
            FOUND_VAR LIBINT2_FOUND
            REQUIRED_VARS LIBINT2_INCLUDE_DIR
            VERSION_VAR LIBINT2_VERSION
            )

    if(LIBINT2_FOUND)
        set(LIBINT2_LIBRARIES ${LIBINT2_LIBRARY})
        set(LIBINT2_INCLUDE_DIRS ${LIBINT2_INCLUDE_DIR} ${LIBINT2_INCLUDE_DIR}/libint2 ${PC_LIBINT2_INCLUDE_DIRS})
        # sanitize LIBINT2_INCLUDE_DIRS: remove duplicates and non-existent entries
        list(REMOVE_DUPLICATES LIBINT2_INCLUDE_DIRS)
        set(LIBINT2_INCLUDE_DIRS_SANITIZED )
        foreach(DIR IN LISTS LIBINT2_INCLUDE_DIRS)
            if (EXISTS ${DIR})
                list(APPEND LIBINT2_INCLUDE_DIRS_SANITIZED ${DIR})
            endif()
        endforeach()
        set(LIBINT2_INCLUDE_DIRS ${LIBINT2_INCLUDE_DIRS_SANITIZED})
    endif()

    if(LIBINT2_FOUND AND NOT TARGET Libint2::Libint)
        add_library(Libint2::int2 INTERFACE IMPORTED)
        set_target_properties(Libint2::int2 PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${LIBINT2_INCLUDE_DIR};${LIBINT2_INCLUDE_DIR}/libint2"
                )
        set_target_properties(Libint2::int2 PROPERTIES
                INTERFACE_LINK_LIBRARIES ${LIBINT2_LIBRARY}
                )
        set_target_properties(Libint2::int2 PROPERTIES
                INTERFACE_COMPILE_FEATURES "cxx_std_11"
                )
    endif()

    if(LIBINT2_FOUND AND NOT TARGET Libint2::cxx)
        add_library(Libint2::cxx INTERFACE IMPORTED)
        set_target_properties(Libint2::cxx PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${LIBINT2_INCLUDE_DIRS}"
                )
        target_link_libraries(Libint2::cxx INTERFACE Libint2::int2)
    endif()

else(PKG_CONFIG_FOUND)

    message(FATAL_ERROR "Could not find the required pkg-config executable")

endif(PKG_CONFIG_FOUND)
