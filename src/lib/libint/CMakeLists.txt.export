cmake_minimum_required(VERSION 3.16)
cmake_policy(SET CMP0074 NEW)

############################# Version and Metadata #############################

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
include(int_computed)

project(
  Libint2
  VERSION ${LIBINT_VERSION}
  DESCRIPTION "High-performance library for computing Gaussian integrals in quantum mechanics"
  HOMEPAGE_URL "https://github.com/evaleev/libint"
  LANGUAGES CXX
  )

set(${PROJECT_NAME}_AUTHORS "Edward F. Valeev")
set(${PROJECT_NAME}_LICENSE "LGPL-3.0")

set(L2 Libint2)  # Namespace
set(pnv libint2) # projectnameversion

################################### Overview ###################################

# * high MAX_AM generating >20k files may require `ulimit -s 65535` for linking
#   library target on Linux to avert "ld: Argument list too long" w/o UNITY_BUILD
# * when building for packaging, absence or presence of Boost determines whether
#   libint2/boost headers are installed or not

################################### Options ####################################
include(GNUInstallDirs)
include(CTest)
include(options)

#  <<<  General  >>>

option_with_default(CMAKE_BUILD_TYPE "Build type" Release)
option_with_print(LIBINT2_REQUIRE_CXX_API "C++11 Libint API: define library targets + test (requires Eigen3)" ON)
option_with_print(LIBINT2_REQUIRE_CXX_API_COMPILED "Build C++11 Compiled (not just header-only) targets (requires Eigen3)" ON)
option_with_print(LIBINT2_ENABLE_FORTRAN "Build Fortran03+ Libint interface (requires C and Fortran and Python)" OFF)
option_with_print(LIBINT2_ENABLE_MPFR "Use GNU MPFR library for high-precision testing (requires MPFR. EXPERTS ONLY)" OFF)
option_with_print(LIBINT2_ENABLE_PYTHON "Build Python bindings (requires Python and Pybind11 and Eigen3)" OFF)
option_with_print(LIBINT2_PREFIX_PYTHON_INSTALL "For LIBINT2_ENABLE_PYTHON=ON, whether to install the Python module in the Linux manner to CMAKE_INSTALL_PREFIX or to not install it. See target libint2-python-wheel for alternate installation in the Python manner to Python_EXECUTABLE's site-packages." OFF)
option_with_print(BUILD_SHARED_LIBS "Build Libint library as shared, not static" OFF)
option_with_print(LIBINT2_BUILD_SHARED_AND_STATIC_LIBS "Build both shared and static Libint libraries in one shot. Uses -fPIC." OFF)
option_with_print(LIBINT2_LOCAL_Eigen3_INSTALL "Install an exported target with hard-coded Eigen3 dependency paths. This is potentially useful and important when consuming the compiled C++11 interface library so that the Libint library build and Libint consumer build use the same Eigen3 installation & ABI. This is at most a convenience when consuming the header-only C++11 interface library. In consumer build, set `LIBINT2_LOCAL_Eigen3_FIND=ON` before `find_package(Libint2) to load the exported Eigen3." OFF)
option_with_default(LIBINT2_REALTYPE "Specifies the floating-point data type used by the library." double)
option_with_print(CMAKE_DISABLE_FIND_PACKAGE_Boost "When Boost required for C++11 API, disable its detection, thereby forcing use of bundled Boost (Standard CMake variable: https://cmake.org/cmake/help/latest/variable/CMAKE_DISABLE_FIND_PACKAGE_PackageName.html)" OFF)
include(xhost)  # defines: option(ENABLE_XHOST "Enable processor-specific optimization" ON)

if ((LIBINT2_REQUIRE_CXX_API_COMPILED OR LIBINT2_ENABLE_PYTHON) AND NOT LIBINT2_REQUIRE_CXX_API)
    set(LIBINT2_REQUIRE_CXX_API 1)
    message(STATUS "Setting option LIBINT2_REQUIRE_CXX_API=ON as needed by LIBINT2_REQUIRE_CXX_API_COMPILED=${LIBINT2_REQUIRE_CXX_API_COMPILED} and/or LIBINT2_ENABLE_PYTHON=${LIBINT2_ENABLE_PYTHON}.")
endif()

if(LIBINT2_ENABLE_FORTRAN)
  include(CheckLanguage)
  check_language(Fortran)
  if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
  else()
    message(FATAL_ERROR "Given LIBINT2_ENABLE_FORTRAN=ON but could not find Fortran compiler. Provide via CMAKE_Fortran_COMPILER")
  endif()
endif()

#  <<<  Ordering Conventions  >>>

# option_with_default(LIBINT2_SHGAUSS_ORDERING "Ordering for shells of solid harmonic Gaussians: standard or gaussian" standard)  # retired Jan 2023
set(LIBINT2_SHGAUSS_ORDERING ${LIBINT2_SHGAUSS_ORDERING})
message(STATUS "Setting option LIBINT2_SHGAUSS_ORDERING: ${LIBINT2_SHGAUSS_ORDERING} (read-only from generation-time)")

set(LIBINT2_CARTGAUSS_ORDERING ${LIBINT2_CARTGAUSS_ORDERING})
message(STATUS "Setting option LIBINT2_CARTGAUSS_ORDERING: ${LIBINT2_CARTGAUSS_ORDERING} (read-only from generation-time)")

set(LIBINT2_SHELL_SET ${LIBINT2_SHELL_SET})
message(STATUS "Setting option LIBINT2_SHELL_SET: ${LIBINT2_SHELL_SET} (read-only from generation-time)")

#  <<<  Miscellaneous  >>>

# next one defined by `include(CTest)`
message(STATUS "Showing option BUILD_TESTING: ${BUILD_TESTING}")

#  <<<  Path  >>>

# * use GNUInstallDirs defaults, but define others for CMake and basis files
# * LIBINT2_ dirs are STRING, not PATH, to preserve relative paths

# next four defined by `include(GNUInstallDirs)`
message(STATUS "Showing option CMAKE_INSTALL_BINDIR: ${CMAKE_INSTALL_BINDIR}")
message(STATUS "Showing option CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "Showing option CMAKE_INSTALL_INCLUDEDIR: ${CMAKE_INSTALL_INCLUDEDIR}")
message(STATUS "Showing option CMAKE_INSTALL_DATADIR: ${CMAKE_INSTALL_DATADIR}")
set(LIBINT2_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${pnv}"
    CACHE STRING "Directory to which CMake files are installed")
message(STATUS "Showing option LIBINT2_INSTALL_CMAKEDIR: ${LIBINT2_INSTALL_CMAKEDIR}")
set(LIBINT2_INSTALL_BASISDIR "${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION}"
    CACHE STRING "Directory to which data files are installed. basis/ directory created within")
message(STATUS "Showing option LIBINT2_INSTALL_BASISDIR: ${LIBINT2_INSTALL_BASISDIR}")
set(LIBINT2_INSTALL_FMODDIR "${CMAKE_INSTALL_INCLUDEDIR}/${pnv}/fortran2/modules"
    CACHE STRING "Directory to which Fortran module files are installed")

######################## Process & Validate Options ###########################
include(autocmake_safeguards)
include(CheckFunctionExists)
include(int_orderings)

check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
if (HAVE_POSIX_MEMALIGN)
    set(LIBINT2_ALIGN_SIZE "0" CACHE STRING "(EXPERT) if posix_memalign is available, this will specify alignment of Libint data, in units of sizeof(LIBINT2_REALTYPE). Default is to use built-in heuristics")
endif()
if (MSVC)
    set(HAVE_POSIX_MEMALIGN 1)  # aliased it
    set(LIBINT2_ALIGN_SIZE "0" CACHE STRING "(LAB) I so don't understand the option handling of ALIGN_SIZE")
endif()

find_package(Eigen3 MODULE)

if (TARGET Eigen3::Eigen)
    set(LIBINT_HAS_EIGEN 1)
    set(LIBINT_HAS_CXX_API 1)
else()
    set(LIBINT_HAS_CXX_API 0)
endif()
if (LIBINT2_REQUIRE_CXX_API AND NOT ${LIBINT_HAS_CXX_API})
    message(FATAL_ERROR "C++ API cannot be built without Eigen3; configure (via CMake) and install Eigen3 and add the install prefix to CMAKE_PREFIX_PATH, or add -DLIBINT2_REQUIRE_CXX_API=OFF to the CMake command line if the C++ API is not required")
endif()

if (LIBINT2_ENABLE_MPFR)
    find_package(Multiprecision MODULE REQUIRED COMPONENTS gmpxx mpfr)

    get_property(_loc TARGET Multiprecision::gmp PROPERTY LOCATION)
    message(STATUS "${Cyan}Found GMP${ColourReset}: ${_loc}")
    get_property(_loc TARGET Multiprecision::gmpxx PROPERTY LOCATION)
    message(STATUS "${Cyan}Found GMPXX${ColourReset}: ${_loc}")
    get_property(_loc TARGET Multiprecision::mpfr PROPERTY LOCATION)
    message(STATUS "${Cyan}Found MPFR${ColourReset}: ${_loc} (found version ${MPFR_VERSION})")
endif()

if (LIBINT2_REQUIRE_CXX_API)
    if (NOT TARGET Boost::headers)
        find_package(Boost 1.57)
    endif()
    if (TARGET Boost::headers)
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

#  <<<  Build  >>>

configure_file(include/libint2/config2.h.cmake.in include/libint2/config2.h @ONLY)
    # configuration.h/cc define a string summary of capabilities. @ONLY is maximally deferred in case config2 changes the
    #   summary (as it used to when LIBINT_SHGSHELL_ORDERING was library-config-time selected).
configure_file(include/libint2/util/configuration.h.cmake.in include/libint2/util/configuration.h @ONLY)
configure_file(src/configuration.cc.cmake.in src/configuration.cc @ONLY)
set(DATADIR_ABSOLUTE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION})
configure_file(include/libint2/basis.h.cmake.in include/libint2/basis.h @ONLY)

include(srclist.cmake)
list(LENGTH LIBINT2_LIBRARY_CXX_SRC _source_count)
message(STATUS "Loading ${_source_count} library source files from LIBINT2_LIBRARY_CXX_SRC")
message(DEBUG "LIBINT2_LIBRARY_CXX_SRC=${LIBINT2_LIBRARY_CXX_SRC}")

set(export_properties
  "Libint2_VERSION"
  "Libint2_MAX_AM_ERI"
  "Libint2_CONFIGURATION"
  )

add_library(int-obj OBJECT "${LIBINT2_LIBRARY_CXX_SRC}" "src/configuration.cc")

target_compile_definitions(int-obj PRIVATE -D__COMPILING_LIBINT2=1)

target_compile_features(int-obj PUBLIC "cxx_std_11") # N.B. PUBLIC to make int-{static/shared} require C++11?

set_target_properties(int-obj PROPERTIES UNITY_BUILD TRUE) # always use unity build for int-obj

if (MSVC)
    # Increase stack size from 1 MB to 4 MB
    set_target_properties(int-obj PROPERTIES LINK_FLAGS "/STACK:4194304")
endif()

target_include_directories(int-obj PRIVATE ${PROJECT_BINARY_DIR}/include/
                                           ${PROJECT_SOURCE_DIR}/src/
                                           ${PROJECT_SOURCE_DIR}/include/
                                           ${PROJECT_SOURCE_DIR}/include/libint2/
                                           )

if(BUILD_SHARED_LIBS OR LIBINT2_BUILD_SHARED_AND_STATIC_LIBS OR (LIBINT2_ENABLE_PYTHON AND NOT MSVC))
    set(BUILD_SHARED_LIBS 1)
    set_target_properties(int-obj PROPERTIES POSITION_INDEPENDENT_CODE 1)
endif()
if (NOT BUILD_SHARED_LIBS OR LIBINT2_BUILD_SHARED_AND_STATIC_LIBS)
    set(BUILD_STATIC_LIBS 1)
endif()

# if building CXX API, make an object lib needed for non-header-only version
if (LIBINT2_REQUIRE_CXX_API_COMPILED)
    add_library(int-cxx-compiled-obj OBJECT src/engine.cpp)

    target_compile_definitions(
      int-cxx-compiled-obj
      PUBLIC
        LIBINT2_DOES_NOT_INLINE_ENGINE=1
        __COMPILING_LIBINT2=1
      )

    target_compile_features(int-cxx-compiled-obj PUBLIC "cxx_std_11")

    if (BUILD_SHARED_LIBS)
        set_target_properties(int-cxx-compiled-obj PROPERTIES POSITION_INDEPENDENT_CODE 1)
    endif()

    if(MSVC)
        # MSVC does not include <cmath> constants, unless _USE_MATH_DEFINES is defined.
        target_compile_definitions(int-cxx-compiled-obj PUBLIC _USE_MATH_DEFINES)
        # Set the exception handling model
        target_compile_options(int-cxx-compiled-obj PUBLIC "/EHsc")
    endif()

    target_include_directories(
      int-cxx-compiled-obj
      PRIVATE
        ${PROJECT_BINARY_DIR}/include/
        ${PROJECT_SOURCE_DIR}/src/
        ${PROJECT_SOURCE_DIR}/include/
        ${PROJECT_BINARY_DIR}/include/libint2/
      )

    target_link_libraries(int-cxx-compiled-obj PRIVATE Eigen3::Eigen)
    if (TARGET Boost::headers)
        target_link_libraries(int-cxx-compiled-obj PRIVATE Boost::headers)
    endif()
endif()

if (BUILD_SHARED_LIBS)
    add_library(int-shared SHARED $<TARGET_OBJECTS:int-obj>)

    target_compile_definitions(
      int-shared
      INTERFACE
        $<BUILD_INTERFACE:__COMPILING_LIBINT2=1>
    )

    if(MSVC)
        target_compile_definitions(int-shared PUBLIC _USE_MATH_DEFINES)
        target_compile_options(int-shared PUBLIC "/EHsc")
    endif()

    target_compile_features(int-shared INTERFACE "cxx_std_11")

    set_target_properties(
      int-shared
      PROPERTIES
        SOVERSION ${LIBINT_MAJOR_SOVERSION}
        MACOSX_RPATH ON
        OUTPUT_NAME "int2"
        EXPORT_NAME "int2"
        Libint2_VERSION "${LIBINT_VERSION}"
        Libint2_MAX_AM_ERI "${Libint2_MAX_AM_ERI}"
        Libint2_CONFIGURATION "${Libint2_CONFIG_COMPONENTS}"
      )
    set_property(TARGET int-shared APPEND PROPERTY EXPORT_PROPERTIES "${export_properties}")
    if (APPLE)
        set_target_properties(
          int-shared
          PROPERTIES
            LINK_FLAGS "-undefined dynamic_lookup"
          )
    endif()
    if (MSVC)
        # necessary but insufficient for dll, https://github.com/evaleev/libint/issues/237
        set_target_properties(
          int-shared
          PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS ON
      )
    endif()

    target_include_directories(
      int-shared
      INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/libint2>
      )

    # C++ library
    if (LIBINT2_REQUIRE_CXX_API)
        add_library(int-cxx-headeronly-shared INTERFACE)

        target_compile_definitions(
          int-cxx-headeronly-shared
          INTERFACE
            $<BUILD_INTERFACE:__COMPILING_LIBINT2=1>
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

        target_compile_features(int-cxx-headeronly-shared INTERFACE "cxx_std_11")

        set_target_properties(int-cxx-headeronly-shared PROPERTIES EXPORT_NAME "cxx")

        target_link_libraries(int-cxx-headeronly-shared INTERFACE int-shared Eigen3::Eigen)
        if (TARGET Boost::headers)
            target_link_libraries(int-cxx-headeronly-shared INTERFACE Boost::headers)
        endif()

        # now make the compiled library
        if (LIBINT2_REQUIRE_CXX_API_COMPILED)
        add_library(int-cxx-compiled-shared SHARED $<TARGET_OBJECTS:int-cxx-compiled-obj>)

        if(MSVC)
            target_compile_definitions(int-cxx-compiled-shared PUBLIC _USE_MATH_DEFINES)
            target_compile_options(int-cxx-compiled-shared PUBLIC "/EHsc")
        endif()

        set_target_properties(
          int-cxx-compiled-shared
          PROPERTIES
            SOVERSION ${LIBINT_MAJOR_SOVERSION}
            MACOSX_RPATH ON
            OUTPUT_NAME "int2-cxx"
            EXPORT_NAME "int2-cxx"
          )
        if (APPLE)
            set_target_properties(
              int-cxx-compiled-shared
              PROPERTIES
                LINK_FLAGS "-undefined dynamic_lookup"
              )
        endif()
        if (MSVC)
            # necessary but insufficient for dll
            set_target_properties(
              int-cxx-compiled-shared
              PROPERTIES
                WINDOWS_EXPORT_ALL_SYMBOLS ON
          )
        endif()

        target_link_libraries(int-cxx-compiled-shared INTERFACE int-cxx-headeronly-shared)
        endif()
    endif()
endif()

if (BUILD_STATIC_LIBS)
    add_library(int-static STATIC $<TARGET_OBJECTS:int-obj>)

    target_compile_features(int-static INTERFACE "cxx_std_11")

    target_compile_definitions(
      int-static
      INTERFACE
        $<BUILD_INTERFACE:__COMPILING_LIBINT2=1>
    )

    if(MSVC)
        target_compile_definitions(int-static PUBLIC _USE_MATH_DEFINES)
        target_compile_options(int-static PUBLIC "/EHsc")
    endif()

    set_target_properties(
      int-static
      PROPERTIES
        OUTPUT_NAME "int2"
        EXPORT_NAME "int2"
        Libint2_VERSION "${LIBINT_VERSION}"
        Libint2_MAX_AM_ERI "${Libint2_MAX_AM_ERI}"
        Libint2_CONFIGURATION "${Libint2_CONFIG_COMPONENTS}"
      )
    set_property(TARGET int-static APPEND PROPERTY EXPORT_PROPERTIES "${export_properties}")

    target_include_directories(
      int-static
      INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/libint2>
      )

    # C++ library
    if (LIBINT2_REQUIRE_CXX_API)
        add_library(int-cxx-headeronly-static INTERFACE)

        target_compile_definitions(
          int-cxx-headeronly-static
          INTERFACE
            $<BUILD_INTERFACE:__COMPILING_LIBINT2=1>
            $<BUILD_INTERFACE:SRCDATADIR="${PROJECT_SOURCE_DIR}/lib/basis">
          )
        if (NOT MSVC)
            target_compile_definitions(
              int-cxx-headeronly-static
              INTERFACE
                $<INSTALL_INTERFACE:DATADIR="\$\{_IMPORT_PREFIX\}/${LIBINT2_INSTALL_BASISDIR}">
              )
        endif()

        target_compile_features(int-cxx-headeronly-static INTERFACE "cxx_std_11")

        set_target_properties(int-cxx-headeronly-static PROPERTIES EXPORT_NAME "cxx")

        target_link_libraries(int-cxx-headeronly-static INTERFACE int-static Eigen3::Eigen)
        if (TARGET Boost::headers)
            target_link_libraries(int-cxx-headeronly-static INTERFACE Boost::headers)
        endif()

        # now make the compiled library
        if (LIBINT2_REQUIRE_CXX_API_COMPILED)
        add_library(int-cxx-compiled-static STATIC $<TARGET_OBJECTS:int-cxx-compiled-obj>)

        if(MSVC)
            target_compile_definitions(int-cxx-compiled-static PUBLIC _USE_MATH_DEFINES)
            target_compile_options(int-cxx-compiled-static PUBLIC "/EHsc")
        endif()

        set_target_properties(
          int-cxx-compiled-static
          PROPERTIES
            OUTPUT_NAME "int2-cxx"
            EXPORT_NAME "int2-cxx"
          )

        target_link_libraries(int-cxx-compiled-static INTERFACE int-cxx-headeronly-static)
        endif()
    endif()
endif()

# Permanent aliases ====================================================================================================
# * make Libint2:: targets available for both `add_subdirectory(Libint2)` & `find_package(Libint2)` approaches
# * used for tests

if (BUILD_SHARED_LIBS)
    add_library(${L2}::int2 ALIAS int-shared)
    if (LIBINT2_REQUIRE_CXX_API)
        add_library(${L2}::cxx ALIAS int-cxx-headeronly-shared)
        if (LIBINT2_REQUIRE_CXX_API_COMPILED)
            add_library(${L2}::int2-cxx ALIAS int-cxx-compiled-shared)
        endif()
    endif()
elseif (BUILD_STATIC_LIBS)
    add_library(${L2}::int2 ALIAS int-static)
    if (LIBINT2_REQUIRE_CXX_API)
        add_library(${L2}::cxx ALIAS int-cxx-headeronly-static)
        if (LIBINT2_REQUIRE_CXX_API_COMPILED)
            add_library(${L2}::int2-cxx ALIAS int-cxx-compiled-static)
        endif()
    endif()
endif()

# Legacy (pre-2.7.0) aliases ===========================================================================================

if (BUILD_SHARED_LIBS)
    add_library(libint2 ALIAS int-shared)
    if (LIBINT2_REQUIRE_CXX_API)
        add_library(libint2_cxx ALIAS int-cxx-headeronly-shared)
    endif()
elseif (BUILD_STATIC_LIBS)
    add_library(libint2 ALIAS int-static)
    if (LIBINT2_REQUIRE_CXX_API)
        add_library(libint2_cxx ALIAS int-cxx-headeronly-static)
    endif()
endif()

# Fortran bindings =====================================================================================================

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
        add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/fortran/libint2.h.i
                COMMAND ${CMAKE_C_COMPILER} -E -D__COMPILING_LIBINT2 -I${PROJECT_SOURCE_DIR}/include -I${PROJECT_SOURCE_DIR}/src -I${PROJECT_BINARY_DIR}/include
                -I${PROJECT_BINARY_DIR}/include/libint2 ${PROJECT_SOURCE_DIR}/include/libint2.h -o ${PROJECT_BINARY_DIR}/fortran/libint2.h.i
                DEPENDS ${PROJECT_SOURCE_DIR}/include/libint2.h
                COMMENT "Generating libint2.h.i"
                )
    else()
        message(FATAL_ERROR "Cannot run preprocessor on non-Unix systems, disable Fortran to proceed")
    endif()

    # translated Libint_t
    add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/fortran/libint2_types_f.h
            COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/fortran/c_to_f.py ${PROJECT_BINARY_DIR}/fortran/libint2.h.i ${PROJECT_BINARY_DIR}/fortran/libint2_types_f.h Libint_t
            DEPENDS ${PROJECT_BINARY_DIR}/fortran/libint2.h.i
            COMMENT "Generating libint2_types_f.h"
            )

    # extracted defines from libint2_types.h
    add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/fortran/fortran_incldefs.h
            COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/fortran/make_defs.py ${PROJECT_SOURCE_DIR}/src/libint2_types.h ${PROJECT_BINARY_DIR}/fortran/fortran_incldefs.h
            DEPENDS ${PROJECT_SOURCE_DIR}/src/libint2_types.h
            COMMENT "Generating fortran_incldefs.h"
            )

    # build module
    add_library(libint_f OBJECT fortran/libint_f.F90)

    set_source_files_properties(fortran/libint_f.F90 PROPERTIES OBJECT_DEPENDS "${PROJECT_BINARY_DIR}/fortran/libint2_types_f.h;${PROJECT_BINARY_DIR}/fortran/fortran_incldefs.h")
    target_compile_definitions(libint_f PRIVATE __COMPILING_LIBINT2)

    set_target_properties(
      libint_f
      PROPERTIES
        Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/${BUILDTREE_FMODDIR}
      )
    if (BUILD_SHARED_LIBS)
        set_target_properties(libint_f PROPERTIES POSITION_INDEPENDENT_CODE 1)
    endif()

    target_include_directories(libint_f PUBLIC
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/fortran>
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${BUILDTREE_FMODDIR}/>
            $<INSTALL_INTERFACE:${LIBINT2_INSTALL_FMODDIR}>)

    # Fortran tests merged into rest of tests

    # install Fortran modules
    install(
      DIRECTORY
        "${PROJECT_BINARY_DIR}/${BUILDTREE_FMODDIR}/"
      COMPONENT "${L2}_Development_Fortran"
      DESTINATION "${LIBINT2_INSTALL_FMODDIR}"
      )

        install(
          TARGETS
            libint_f
          EXPORT fshared_set
          RUNTIME
            COMPONENT ${L2}_Runtime
          LIBRARY
            COMPONENT ${L2}_Runtime
            NAMELINK_COMPONENT ${L2}_Development
          ARCHIVE
            COMPONENT ${L2}_Development
          PUBLIC_HEADER
            COMPONENT ${L2}_Development
          INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
          )

endif()

# Python bindings ======================================================================================================

if (LIBINT2_ENABLE_PYTHON)
    add_subdirectory(python)
endif()

#  <<<  Install  >>>

install(
  DIRECTORY
    ${PROJECT_BINARY_DIR}/include/
    ${PROJECT_SOURCE_DIR}/include/
  TYPE INCLUDE
  COMPONENT ${L2}_Development
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
  )

install(
  FILES
    ${PROJECT_SOURCE_DIR}/src/libint2_iface.h
    ${PROJECT_SOURCE_DIR}/src/libint2_iface_internal.h
    ${PROJECT_SOURCE_DIR}/src/libint2_params.h
    ${PROJECT_SOURCE_DIR}/src/libint2_types.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libint2
  COMPONENT ${L2}_Development
  )

if (BUILD_SHARED_LIBS)
    if (TARGET int-cxx-compiled-shared)
        set(tgts int-shared int-cxx-headeronly-shared int-cxx-compiled-shared)
    elseif (TARGET int-cxx-headeronly-shared)
        set(tgts int-shared int-cxx-headeronly-shared)
    else()
        set(tgts int-shared)
    endif()

    install(
      TARGETS
        ${tgts}
      EXPORT shared_set
      RUNTIME
        COMPONENT ${L2}_Runtime
      LIBRARY
        COMPONENT ${L2}_Runtime
        NAMELINK_COMPONENT ${L2}_Development
      ARCHIVE
        COMPONENT ${L2}_Development
      PUBLIC_HEADER
        COMPONENT ${L2}_Development
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
      )
endif()

if (BUILD_STATIC_LIBS)
    if (TARGET int-cxx-compiled-static)
        set(tgts int-static int-cxx-headeronly-static int-cxx-compiled-static)
    elseif (TARGET int-cxx-headeronly-static)
        set(tgts int-static int-cxx-headeronly-static)
    else()
        set(tgts int-static)
    endif()

    install(
      TARGETS
        ${tgts}
      EXPORT static_set
      RUNTIME
        COMPONENT ${L2}_Runtime
      LIBRARY
        COMPONENT ${L2}_Runtime
        NAMELINK_COMPONENT ${L2}_Development
      ARCHIVE
        COMPONENT ${L2}_Development
      PUBLIC_HEADER
        COMPONENT ${L2}_Development
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
      )
endif()

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

# install basis set library
install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/lib/basis
  DESTINATION "${LIBINT2_INSTALL_BASISDIR}"
  COMPONENT ${L2}_Development
  )

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()

#  <<<  Export Config  >>>

configure_package_config_file(
  cmake/${pnv}-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${pnv}-config.cmake"
  INSTALL_DESTINATION ${LIBINT2_INSTALL_CMAKEDIR}
  NO_SET_AND_CHECK_MACRO
  )
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${pnv}-config-version.cmake
                                 VERSION ${LIBINT_VERSION}
                                 COMPATIBILITY SameMajorVersion)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${pnv}-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${pnv}-config-version.cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/FindEigen3.cmake
  DESTINATION ${LIBINT2_INSTALL_CMAKEDIR}
  COMPONENT ${L2}_Development
  )

include(JoinPaths)
join_paths(libdir_for_pc_file "\${exec_prefix}" "${CMAKE_INSTALL_LIBDIR}")
join_paths(includedir_for_pc_file "\${prefix}" "${CMAKE_INSTALL_INCLUDEDIR}")

configure_file(cmake/libint2.pc.cmake.in libint2.pc @ONLY)
if (NOT MSVC)
    install(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/libint2.pc
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig/
      COMPONENT ${L2}_Development
      )
endif()

if (BUILD_SHARED_LIBS)
    install(
      EXPORT shared_set
      FILE "${pnv}-targets-shared.cmake"
      NAMESPACE "${L2}::"
      DESTINATION ${LIBINT2_INSTALL_CMAKEDIR}
      COMPONENT ${L2}_Development
      )
endif()

if (BUILD_STATIC_LIBS)
    install(
      EXPORT static_set
      FILE "${pnv}-targets-static.cmake"
      NAMESPACE "${L2}::"
      DESTINATION ${LIBINT2_INSTALL_CMAKEDIR}
      COMPONENT ${L2}_Development
      )
endif()

if (LIBINT2_LOCAL_Eigen3_INSTALL AND TARGET Eigen3::Eigen)
    install(
      EXPORT local_set
      FILE "${pnv}-targets-eigen3.cmake"
      NAMESPACE "${L2}::"
      DESTINATION ${LIBINT2_INSTALL_CMAKEDIR}
      COMPONENT ${L2}_Eigen3
      )
endif()

# not much point in exporting Fortran since it's a simple int2+includes+flags and
#   those includes are compiler-version dependent, but this is the setup.
#   note that FILE something like "fortran-shared" NOT "shared-fortran" or the glob
#   confused with shared_set
#
#    if(LIBINT2_ENABLE_FORTRAN)
#        install(
#          EXPORT fshared_set
#          FILE "${pnv}-targets-fortran-shared.cmake"
#          NAMESPACE "${L2}::"
#          DESTINATION ${LIBINT2_INSTALL_CMAKEDIR}
#          COMPONENT ${L2}_Development_Fortran
#          )
#    endif()
