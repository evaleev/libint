cmake_minimum_required(VERSION 3.16)  # 3.16: unity build
                                      # 3.15: new Python detection
                                      # 3.8: introduced C++ standards as features
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0077 NEW)
cmake_policy(SET CMP0079 NEW)

############################# Version and Metadata #############################

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
include(int_computed)

project(
  Libint2
  VERSION ${LIBINT_VERSION}
  DESCRIPTION ${LIBINT_DESCRIPTION}
  HOMEPAGE_URL "http://libint.valeyev.net"
  LANGUAGES CXX
  )

set(${PROJECT_NAME}_AUTHORS "Edward F. Valeev")
set(${PROJECT_NAME}_LICENSE "LGPL-3.0-only")

set(L2 Libint2)  # Namespace
set(pnv libint2) # projectnameversion

#################################### Guide #####################################

# See "TARBALL" labels in INSTALL.md for elaboration of options, dependencies, & targets.

################################### Options ####################################
include(GNUInstallDirs)
include(CTest)
include(options)
include(CheckFunctionExists)

#  <<<  General  >>>

option_with_default(CMAKE_BUILD_TYPE
  "Build type" Release)
option_with_print(BUILD_SHARED_LIBS
  "Build Libint library as shared, not static" OFF)
option_with_print(LIBINT2_BUILD_SHARED_AND_STATIC_LIBS
  "Build both shared and static Libint libraries in one shot. Uses -fPIC." OFF)
option_with_print(LIBINT2_REQUIRE_CXX_API
  "C++11 Libint API: define library targets + test (requires Eigen3, Boost is optional but strongly recommended)" ON)
option_with_print(LIBINT2_REQUIRE_CXX_API_COMPILED
  "Build C++11 Compiled (not just header-only) targets (requires Eigen3)" OFF)
option_with_print(LIBINT2_ENABLE_FORTRAN
  "Build Fortran03+ Libint interface (requires C and Fortran and Python)" OFF)
option_with_print(LIBINT2_ENABLE_MPFR
  "Use GNU MPFR library for high-precision testing (requires MPFR. EXPERTS ONLY)" OFF)
option_with_print(LIBINT2_ENABLE_PYTHON
  "Build Python bindings (requires Python and Pybind11 and Eigen3)" OFF)
option_with_print(LIBINT2_PREFIX_PYTHON_INSTALL
  "For LIBINT2_ENABLE_PYTHON=ON, whether to install the Python module in the Linux manner to CMAKE_INSTALL_PREFIX or to not install it. See target libint2-python-wheel for alternate installation in the Python manner to Python_EXECUTABLE's site-packages." OFF)
option_with_print(LIBINT2_LOCAL_Eigen3_INSTALL
  "Install an exported target with hard-coded Eigen3 dependency paths. This is potentially useful and important when consuming the compiled C++11 interface library so that the Libint library build and Libint consumer build use the same Eigen3 installation & ABI. This is at most a convenience when consuming the header-only C++11 interface library. In consumer build, set `LIBINT2_LOCAL_Eigen3_FIND=ON` before `find_package(Libint2) to load the exported Eigen3." OFF)
option_with_print(CMAKE_DISABLE_FIND_PACKAGE_Boost
  "When Boost required for C++11 API, disable its detection, thereby forcing use of bundled Boost (Standard CMake variable: https://cmake.org/cmake/help/latest/variable/CMAKE_DISABLE_FIND_PACKAGE_PackageName.html)" OFF)

check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
if (HAVE_POSIX_MEMALIGN)
    option_with_default(LIBINT2_ALIGN_SIZE
      "(EXPERT) if posix_memalign is available, this will specify alignment of Libint data, in units of sizeof(LIBINT2_REALTYPE). Default is to use built-in heuristics." 0)
    mark_as_advanced(LIBINT2_ALIGN_SIZE)
elseif (MSVC)
    # works for Windows in memory.h
    set(LIBINT2_ALIGN_SIZE "0" CACHE STRING "")
endif()

option_with_default(LIBINT2_REALTYPE
  "Specifies the floating-point data type used by the library." double)
option_with_print(LIBINT_USER_DEFINED_REAL_INCLUDES
  "UNTESTED Additional #includes necessary to use the real type." OFF)
include(int_userreal)

if ((LIBINT2_REQUIRE_CXX_API_COMPILED OR LIBINT2_ENABLE_PYTHON) AND NOT LIBINT2_REQUIRE_CXX_API)
    set(LIBINT2_REQUIRE_CXX_API 1)
    message(STATUS "Setting option LIBINT2_REQUIRE_CXX_API=ON as needed by LIBINT2_REQUIRE_CXX_API_COMPILED=${LIBINT2_REQUIRE_CXX_API_COMPILED} and/or LIBINT2_ENABLE_PYTHON=${LIBINT2_ENABLE_PYTHON}.")
endif()

if (LIBINT2_ENABLE_FORTRAN)
    include(CheckLanguage)
    check_language(Fortran)
    if (CMAKE_Fortran_COMPILER)
        enable_language(Fortran)
    else()
        message(FATAL_ERROR "Given LIBINT2_ENABLE_FORTRAN=ON but could not find Fortran compiler. Provide via CMAKE_Fortran_COMPILER")
    endif()
endif()

#  <<<  Miscellaneous  >>>

# next one defined by `include(CTest)`
message(STATUS "Showing option BUILD_TESTING: ${BUILD_TESTING}")

#  <<<  Ordering Conventions  >>>

# Jan 2023, retired  # option_with_default(LIBINT2_SHGAUSS_ORDERING ...
# August 2023: the following variable has an effect on `INT_SOLIDHARMINDEX(l, m)` but otherwise the choice can be deferred to runtime.
set(LIBINT2_SHGAUSS_ORDERING ${LIBINT2_SHGAUSS_ORDERING})
message(STATUS "Setting option LIBINT2_SHGAUSS_ORDERING: ${LIBINT2_SHGAUSS_ORDERING} (read-only from generation-time)")

set(LIBINT2_CARTGAUSS_ORDERING ${LIBINT2_CARTGAUSS_ORDERING})
message(STATUS "Setting option LIBINT2_CARTGAUSS_ORDERING: ${LIBINT2_CARTGAUSS_ORDERING} (read-only from generation-time)")

set(LIBINT2_SHELL_SET ${LIBINT2_SHELL_SET})
message(STATUS "Setting option LIBINT2_SHELL_SET: ${LIBINT2_SHELL_SET} (read-only from generation-time)")


######################## Process & Validate Options ###########################
include(autocmake_safeguards)


# ====  Dependencies  ===========================================================

# See notes at https://github.com/evaleev/libint/blob/master/INSTALL.md#prerequisites

if (LIBINT2_ENABLE_MPFR)
    find_package(Multiprecision MODULE REQUIRED COMPONENTS gmpxx mpfr)
    set(LIBINT_HAS_MPFR 1)

    get_property(_loc TARGET Multiprecision::gmp PROPERTY LOCATION)
    message(VERBOSE "${Cyan}Found GMP${ColourReset}: ${_loc}")
    get_property(_loc TARGET Multiprecision::gmpxx PROPERTY LOCATION)
    message(VERBOSE "${Cyan}Found GMPXX${ColourReset}: ${_loc}")
    get_property(_loc TARGET Multiprecision::mpfr PROPERTY LOCATION)
    message(VERBOSE "${Cyan}Found MPFR${ColourReset}: ${_loc} (found version ${MPFR_VERSION})")
endif()

if (LIBINT2_REQUIRE_CXX_API)
    if (NOT TARGET Boost::headers)
        find_package(Boost 1.57)
    endif()
    if (TARGET Boost::headers)
        include(int_checkboost)
        set(LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS 1)
    else()
        set(LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS 0)

        file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/include/libint2)
        execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/external/boost.tar.gz
          WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/include/libint2
          RESULT_VARIABLE UNPACK_BOOST_RESULT
          OUTPUT_VARIABLE UNPACK_BOOST_OUTPUT
          ERROR_VARIABLE UNPACK_BOOST_OUTPUT
          )
        message(STATUS "Unpacking bundled Boost")
        if (NOT UNPACK_BOOST_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to unpack the bundled Boost! The tar command output:\n${UNPACK_BOOST_OUTPUT}")
        endif()
    endif()
endif()

find_package(Eigen3 MODULE)

if (TARGET Eigen3::Eigen)
    set(LIBINT_HAS_EIGEN 1)
endif()
if (LIBINT2_REQUIRE_CXX_API AND NOT ${LIBINT_HAS_EIGEN})
    message(FATAL_ERROR "C++ API cannot be built without Eigen3; configure (via CMake) and install Eigen3 and add the install prefix to CMAKE_PREFIX_PATH, or add -D LIBINT2_REQUIRE_CXX_API=OFF to the CMake command line if the C++ API is not required")
endif()

# Python is optionally used for testing.
# * But for Fortran, it's additionally required for preprocessing.
# * And for Python bindings, it's required along with its headers.
if (LIBINT2_ENABLE_PYTHON)
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
elseif (LIBINT2_ENABLE_FORTRAN)
    find_package(Python COMPONENTS Interpreter REQUIRED)
else()
    find_package(Python COMPONENTS Interpreter)
endif()


################################# Main Project #################################
include(CMakePackageConfigHelpers)

configure_file(
  include/libint2/config2.h.cmake.in
  include/libint2/config2.h
  @ONLY)
configure_file(
  # configuration.cc defines a string summary of capabilities. @ONLY is maximally
  #   deferred in case config2 changes the summary (as it used to when
  #   LIBINT_SHGSHELL_ORDERING was library-config-time selected).
  src/configuration.cc.cmake.in
  src/configuration.cc
  @ONLY)

set(DATADIR_ABSOLUTE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION})
configure_file(
  include/libint2/basis.h.in  # temp until libtool dropped and moved to basis.h.cmake.in
  include/libint2/basis.h
  @ONLY)

include(srclist.cmake)
list(LENGTH LIBINT2_LIBRARY_CXX_SRC _source_count)
message(STATUS "Loading ${_source_count} library source files from LIBINT2_LIBRARY_CXX_SRC")
message(DEBUG "LIBINT2_LIBRARY_CXX_SRC=${LIBINT2_LIBRARY_CXX_SRC}")

# Notes:
# * __COMPILING_LIBINT2 is needed whenever using libint targets in the build tree
# * MSVC does not include <cmath> constants, unless _USE_MATH_DEFINES is defined
# * /EHsc sets the exception handling model (allows "throw")
# * _CRT_* squashs some getenv, strdup, strncpy, ctime, fopen warnings

# ====  pre-user-targets: plain/C++(compiled)  ==================================

add_library(
  int-obj
  OBJECT
    "${LIBINT2_LIBRARY_CXX_SRC}"
    src/configuration.cc
  )
target_compile_features(
  int-obj
  PUBLIC
    "cxx_std_11"  # N.B. PUBLIC to make int-{static/shared} require C++11
  )
target_compile_definitions(
  int-obj
  PRIVATE
    __COMPILING_LIBINT2=1
  )
set_target_properties(
  int-obj
  PROPERTIES
    UNITY_BUILD TRUE  # always use unity build for int-obj
  )
if (MSVC)
    set_target_properties(
      int-obj
      PROPERTIES
        LINK_FLAGS "/STACK:4194304"
        # Increase stack size from 1 MB to 4 MB
      )
endif()
target_include_directories(
  int-obj
  PRIVATE
    ${PROJECT_SOURCE_DIR}/include/
    ${PROJECT_BINARY_DIR}/include/
  )

if (LIBINT2_REQUIRE_CXX_API_COMPILED)

    add_library(
      int-cxx-obj
      OBJECT
        src/engine.cpp
      )
    target_compile_features(
      int-cxx-obj
      PUBLIC
        "cxx_std_11"
      )
    target_compile_definitions(
      int-cxx-obj
      PRIVATE
        __COMPILING_LIBINT2=1
        LIBINT2_DOES_NOT_INLINE_ENGINE=1
        $<$<BOOL:${MSVC}>:_USE_MATH_DEFINES>
        $<$<BOOL:${MSVC}>:_CRT_NONSTDC_NO_DEPRECATE>
        $<$<BOOL:${MSVC}>:_CRT_NONSTDC_NO_WARNINGS>
        $<$<BOOL:${MSVC}>:_CRT_SECURE_NO_WARNINGS>
      )
    target_include_directories(
      int-cxx-obj
      PRIVATE
        ${PROJECT_SOURCE_DIR}/include/
        ${PROJECT_BINARY_DIR}/include/
      )
    target_link_libraries (
      int-cxx-obj
      PRIVATE
        Eigen3::Eigen
        $<$<TARGET_EXISTS:Boost::headers>:Boost::headers>
      )

endif()

# plan shared, static, or both, then set fpic accordingly

if (LIBINT2_BUILD_SHARED_AND_STATIC_LIBS OR BUILD_SHARED_LIBS)
    set(L2_BUILD_SHARED_LIBS 1)
endif()
if (LIBINT2_BUILD_SHARED_AND_STATIC_LIBS OR (NOT BUILD_SHARED_LIBS))
    set(L2_BUILD_STATIC_LIBS 1)
endif()
if (L2_BUILD_SHARED_LIBS OR (LIBINT2_ENABLE_PYTHON AND NOT MSVC))
    if (TARGET int-cxx-obj)
        set(tgts int-obj int-cxx-obj)
    else()
        set(tgts int-obj)
    endif()

    set_target_properties(
      ${tgts}
      PROPERTIES
        POSITION_INDEPENDENT_CODE 1
      )
endif()


# ====  6 user targets: plain/C++(headers)/C++(compiled) shared/static  =========

if (L2_BUILD_SHARED_LIBS)

    add_library(
      int-shared
      SHARED
        $<TARGET_OBJECTS:int-obj>
      )
    target_compile_features(
      int-shared
      INTERFACE
        "cxx_std_11"
      )
    target_compile_definitions(
      int-shared
      INTERFACE
        $<BUILD_INTERFACE:__COMPILING_LIBINT2=1>
      PUBLIC
        $<$<BOOL:${MSVC}>:_USE_MATH_DEFINES>
        $<$<BOOL:${MSVC}>:_CRT_NONSTDC_NO_DEPRECATE>
        $<$<BOOL:${MSVC}>:_CRT_NONSTDC_NO_WARNINGS>
        $<$<BOOL:${MSVC}>:_CRT_SECURE_NO_WARNINGS>
      )
    target_compile_options(
      int-shared
      PUBLIC
        $<$<BOOL:${MSVC}>:/EHsc>
      )
    set_target_properties(
      int-shared
      PROPERTIES
        OUTPUT_NAME "int2"
      )
    if (APPLE)
        set_target_properties(
          int-shared
          PROPERTIES
            LINK_FLAGS "-undefined dynamic_lookup"
          )
    endif()
    target_include_directories(
      int-shared
      INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      )

    if (LIBINT2_REQUIRE_CXX_API)

        add_library(
          int-cxx-headeronly-shared
          INTERFACE
          )
        target_compile_features(
          int-cxx-headeronly-shared
          INTERFACE
            "cxx_std_11"
          )
        target_compile_definitions(
          int-cxx-headeronly-shared
          INTERFACE
            $<BUILD_INTERFACE:SRCDATADIR="${PROJECT_SOURCE_DIR}/lib/basis">
          )
        if (NOT MSVC)
            # TODO fix the DATADIR define escaping on Windows
            # * below works fine in tests
            # * but fails in Psi4 compile
            # * prefix replacement in conda used instead on Windows
            target_compile_definitions(
              int-cxx-headeronly-shared
              INTERFACE
                $<INSTALL_INTERFACE:DATADIR="\$\{_IMPORT_PREFIX\}/${LIBINT2_INSTALL_BASISDIR}">
              )
        endif()
        target_link_libraries(
          int-cxx-headeronly-shared
          INTERFACE
            int-shared
            Eigen3::Eigen
            $<$<TARGET_EXISTS:Boost::headers>:Boost::headers>
          )

        if (LIBINT2_REQUIRE_CXX_API_COMPILED)

            add_library(
              int-cxx-shared
              SHARED
                $<TARGET_OBJECTS:int-cxx-obj>
              )
            target_compile_definitions(
              int-cxx-shared
              INTERFACE
                $<BUILD_INTERFACE:__COMPILING_LIBINT2=1>
              PUBLIC
                LIBINT2_DOES_NOT_INLINE_ENGINE=1
                LIBINT2_CONSTEXPR_STATICS=1  # LAB: needed for tests, but correct in general, EFV?
                $<$<BOOL:${MSVC}>:_USE_MATH_DEFINES>
              )
            target_compile_options(
              int-cxx-shared
              PUBLIC
                $<$<BOOL:${MSVC}>:/EHsc>
              )
            set_target_properties(
              int-cxx-shared
              PROPERTIES
                #SOVERSION ${LIBINT_MAJOR_SOVERSION}
                MACOSX_RPATH ON
                OUTPUT_NAME "int2-cxx"
              )
            if (APPLE)
                set_target_properties(
                  int-cxx-shared
                  PROPERTIES
                    LINK_FLAGS "-undefined dynamic_lookup"
                  )
            endif()
            target_link_libraries(
              int-cxx-shared
              INTERFACE
                int-cxx-headeronly-shared
              )

        endif (LIBINT2_REQUIRE_CXX_API_COMPILED)
    endif (LIBINT2_REQUIRE_CXX_API)
endif (L2_BUILD_SHARED_LIBS)

if (L2_BUILD_STATIC_LIBS)

    add_library(
      int-static
      STATIC
        $<TARGET_OBJECTS:int-obj>
      )
    target_compile_features(
      int-static
      INTERFACE
        "cxx_std_11"
      )
    target_compile_definitions(
      int-static
      INTERFACE
        $<BUILD_INTERFACE:__COMPILING_LIBINT2=1>
      PUBLIC
        $<$<BOOL:${MSVC}>:_USE_MATH_DEFINES>
        $<$<BOOL:${MSVC}>:_CRT_NONSTDC_NO_DEPRECATE>
        $<$<BOOL:${MSVC}>:_CRT_NONSTDC_NO_WARNINGS>
        $<$<BOOL:${MSVC}>:_CRT_SECURE_NO_WARNINGS>
      )
    target_compile_options(
      int-static
      PUBLIC
        $<$<BOOL:${MSVC}>:/EHsc>
      )
    set_target_properties(
      int-static
      PROPERTIES
        OUTPUT_NAME "int2"
      )
    target_include_directories(
      int-static
      INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      )

    if (LIBINT2_REQUIRE_CXX_API)

        add_library(
          int-cxx-headeronly-static
          INTERFACE
          )
        target_compile_features(
          int-cxx-headeronly-static
          INTERFACE
            "cxx_std_11"
          )
        target_compile_definitions(
          int-cxx-headeronly-static
          INTERFACE
            $<BUILD_INTERFACE:SRCDATADIR="${PROJECT_SOURCE_DIR}/lib/basis">
          )
        target_link_libraries(
          int-cxx-headeronly-static
          INTERFACE
            int-static
            Eigen3::Eigen
            $<$<TARGET_EXISTS:Boost::headers>:Boost::headers>
          )

        if (LIBINT2_REQUIRE_CXX_API_COMPILED)

            add_library(
              int-cxx-static
              STATIC
                $<TARGET_OBJECTS:int-cxx-obj>
              )
            target_compile_definitions(
              int-cxx-static
              INTERFACE
                $<BUILD_INTERFACE:__COMPILING_LIBINT2=1>
              PUBLIC
                LIBINT2_DOES_NOT_INLINE_ENGINE=1
                LIBINT2_CONSTEXPR_STATICS=1  # LAB: needed for tests, but correct in general, EFV?
                $<$<BOOL:${MSVC}>:_USE_MATH_DEFINES>
              )
            target_compile_options(
              int-cxx-static
              PUBLIC
                $<$<BOOL:${MSVC}>:/EHsc>
              )
            set_target_properties(
              int-cxx-static
              PROPERTIES
                OUTPUT_NAME "int2-cxx"
                #EXPORT_NAME "int2-cxx"
              )
            target_link_libraries(
              int-cxx-static
              INTERFACE
                int-cxx-headeronly-static
              )

        endif (LIBINT2_REQUIRE_CXX_API_COMPILED)
    endif (LIBINT2_REQUIRE_CXX_API)
endif (L2_BUILD_STATIC_LIBS)


# ====  aliases  ================================================================

# permanent aliases
# * used for tests

if (L2_BUILD_SHARED_LIBS)
    add_library(${L2}::int2 ALIAS int-shared)
    if (LIBINT2_REQUIRE_CXX_API)
        add_library(${L2}::cxx ALIAS int-cxx-headeronly-shared)
        if (LIBINT2_REQUIRE_CXX_API_COMPILED)
            add_library(${L2}::int2-cxx ALIAS int-cxx-shared)
        endif()
    endif()
elseif (L2_BUILD_STATIC_LIBS)
    add_library(${L2}::int2 ALIAS int-static)
    if (LIBINT2_REQUIRE_CXX_API)
        add_library(${L2}::cxx ALIAS int-cxx-headeronly-static)
        if (LIBINT2_REQUIRE_CXX_API_COMPILED)
            add_library(${L2}::int2-cxx ALIAS int-cxx-static)
        endif()
    endif()
endif()

# legacy (pre-2.9.0) aliases

if (L2_BUILD_SHARED_LIBS)
    add_library(libint2 ALIAS int-shared)
    if (LIBINT2_REQUIRE_CXX_API)
        add_library(libint2_cxx ALIAS int-cxx-headeronly-shared)
    endif()
elseif (L2_BUILD_STATIC_LIBS)
    add_library(libint2 ALIAS int-static)
    if (LIBINT2_REQUIRE_CXX_API)
        add_library(libint2_cxx ALIAS int-cxx-headeronly-static)
    endif()
endif()


# ====  Fortran bindings  =======================================================

if (LIBINT2_ENABLE_FORTRAN)
    # specify the location of modules
    set(BUILDTREE_FMODDIR "fortran/modules")

    # preprocess libint2.h ... this is a guess for UNIX systems only
    # N.B. Requires C compiler!
    if (UNIX)
        file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/fortran)

        include(CheckLanguage)
        check_language(C)
        if (CMAKE_C_COMPILER)
            enable_language(C)
        else()
            message(FATAL_ERROR "Given LIBINT2_ENABLE_FORTRAN=ON but could not find C compiler needed to generate Fortran bindings, provide via CMAKE_C_COMPILER")
        endif()

        # preprocessed libint.h
        add_custom_command(
          OUTPUT
            ${PROJECT_BINARY_DIR}/fortran/libint2.h.i
          COMMAND ${CMAKE_C_COMPILER} -E
                  -D__COMPILING_LIBINT2
                  -I${PROJECT_SOURCE_DIR}/include
                  -I${PROJECT_SOURCE_DIR}/src
                  -I${PROJECT_BINARY_DIR}/include
                  -I${PROJECT_BINARY_DIR}/include/libint2
                  ${PROJECT_SOURCE_DIR}/include/libint2.h
                  -o ${PROJECT_BINARY_DIR}/fortran/libint2.h.i
          DEPENDS
            ${PROJECT_SOURCE_DIR}/include/libint2.h
          COMMENT "Generating libint2.h.i"
          )
    else()
        message(FATAL_ERROR "Cannot run preprocessor on non-Unix systems, disable Fortran to proceed")
    endif()

    # translated Libint_t
    add_custom_command(
      OUTPUT
        ${PROJECT_BINARY_DIR}/fortran/libint2_types_f.h
      COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/fortran/c_to_f.py ${PROJECT_BINARY_DIR}/fortran/libint2.h.i ${PROJECT_BINARY_DIR}/fortran/libint2_types_f.h Libint_t
      DEPENDS
        ${PROJECT_BINARY_DIR}/fortran/libint2.h.i
      COMMENT "Generating libint2_types_f.h"
      )

    # extracted defines from libint2_types.h
    add_custom_command(
      OUTPUT
        ${PROJECT_BINARY_DIR}/fortran/fortran_incldefs.h
      COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/fortran/make_defs.py ${PROJECT_SOURCE_DIR}/include/libint2_types.h ${PROJECT_BINARY_DIR}/fortran/fortran_incldefs.h
      DEPENDS
        ${PROJECT_SOURCE_DIR}/include/libint2_types.h
      COMMENT "Generating fortran_incldefs.h"
      )

    # build module
    add_library(
      libint_f
      OBJECT
        fortran/libint_f.F90
      )
    set_source_files_properties(
      fortran/libint_f.F90
      PROPERTIES
        OBJECT_DEPENDS "${PROJECT_BINARY_DIR}/fortran/libint2_types_f.h;${PROJECT_BINARY_DIR}/fortran/fortran_incldefs.h"
      )
    target_compile_definitions(
      libint_f
      PRIVATE
        __COMPILING_LIBINT2
      )
    set_target_properties(
      libint_f
      PROPERTIES
        Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/${BUILDTREE_FMODDIR}
      )
    if (BUILD_SHARED_LIBS)
        set_target_properties(libint_f PROPERTIES POSITION_INDEPENDENT_CODE 1)
    endif()
    target_include_directories(
      libint_f
      PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/fortran>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${BUILDTREE_FMODDIR}/>
      )

    # Fortran tests merged into rest of tests

endif (LIBINT2_ENABLE_FORTRAN)


# ====  Python bindings  ========================================================

if (LIBINT2_ENABLE_PYTHON)
    add_subdirectory(python EXCLUDE_FROM_ALL)
endif()


if (BUILD_TESTING)
    add_subdirectory(tests)
endif()

configure_package_config_file(
  cmake/${pnv}-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${pnv}-config.cmake"
  INSTALL_DESTINATION "lib/cmake/libint2" # TODO ${LIBINT2_INSTALL_CMAKEDIR}
  NO_SET_AND_CHECK_MACRO
  )

#  <<<  Install  >>>

install(
  DIRECTORY
    ${PROJECT_SOURCE_DIR}/include2/  # TODO temp for skeleton
  TYPE INCLUDE
  COMPONENT ${L2}_Development
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
  )

if (LIBINT2_LOCAL_Eigen3_INSTALL AND TARGET Eigen3::Eigen)

    add_library(Eigen INTERFACE)
    foreach(prop
      INTERFACE_INCLUDE_DIRECTORIES
      INTERFACE_COMPILE_DEFINITIONS
      INTERFACE_COMPILE_OPTIONS
      INTERFACE_LINK_LIBRARIES
      INTERFACE_POSITION_INDEPENDENT_CODE
      )
        get_property(_propval TARGET Eigen3::Eigen PROPERTY ${prop})
        set_property(TARGET Eigen PROPERTY ${prop} ${_propval})
    endforeach()

    install(
      TARGETS
        Eigen
      EXPORT local_set
      LIBRARY
        COMPONENT ${L2}_Eigen3
        NAMELINK_COMPONENT ${L2}_Eigen3
      ARCHIVE
        COMPONENT ${L2}_Eigen3
      PUBLIC_HEADER
        COMPONENT ${L2}_Eigen3
      )
endif()

#  <<<  Export Config  >>>

if (LIBINT2_LOCAL_Eigen3_INSTALL AND TARGET Eigen3::Eigen)
    install(
      EXPORT local_set
      FILE "${pnv}-targets-eigen3.cmake"
      NAMESPACE "${L2}::"
      DESTINATION ${LIBINT2_INSTALL_CMAKEDIR}
      COMPONENT ${L2}_Eigen3
      )
endif()
