cmake_minimum_required(VERSION 3.8) # introduced C++ standards as features

project(Libint LANGUAGES CXX C)

# Set Libint version ===================================================================================================
set(LIBINT_MAJOR_VERSION 2)
set(LIBINT_MINOR_VERSION 7)
set(LIBINT_MICRO_VERSION 0)
set(LIBINT_BUILDID "beta.3")
set(LIBINT_VERSION "${LIBINT_MAJOR_VERSION}.${LIBINT_MINOR_VERSION}.${LIBINT_MICRO_VERSION}")
if (LIBINT_BUILDID)
    set(LIBINT_EXT_VERSION "${LIBINT_MAJOR_VERSION}.${LIBINT_MINOR_VERSION}.${LIBINT_MICRO_VERSION}-${LIBINT_BUILDID}")
else(LIBINT_BUILDID)
    set(LIBINT_EXT_VERSION ${LIBINT_VERSION})
endif(LIBINT_BUILDID)

# Add module directory and modules =====================================================================================

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules/")
include(CheckCXXSourceCompiles)
include(CheckFunctionExists)
include(CMakeDependentOption)
include(CMakePackageConfigHelpers)
include(FeatureSummary)
include(RedefaultableOption)
include(CMakePushCheckState)

# Options ==============================================================================================================

redefaultable_option(REQUIRE_CXX_API "Fortran03+ Libint interface" ON)
redefaultable_option(ENABLE_FORTRAN "Fortran03+ Libint interface" OFF)
redefaultable_option(ENABLE_MPFR "Use GNU MPFR library for high-precision testing (EXPERTS ONLY)" OFF)
redefaultable_option(LIBINT2_BUILD_SHARED_AND_STATIC_LIBS "Build in one shot static and shared library variants. Uses -fPIC." OFF)
if (NOT LIBINT2_REALTYPE)
    set(LIBINT2_REALTYPE double)
endif()

if (ENABLE_FORTRAN)
    enable_language(Fortran)
endif(ENABLE_FORTRAN)

# SHG ordering
# these are known orderings, must match config.h
set(LIBINT_SHGSHELL_ORDERING_STANDARD 1)
set(LIBINT_SHGSHELL_ORDERING_GAUSSIAN 2)
set(LIBINT2_SHGAUSS_ORDERING "standard" CACHE STRING "Use one of the following known orderings for shells of solid harmonic Gaussians:
  standard -- standard ordering (-l, -l+1 ... l)
  gaussian -- the Gaussian ordering (0, 1, -1, 2, -2, ... l, -l)
The default is standard.")
if (LIBINT2_SHGAUSS_ORDERING STREQUAL "standard")
    set(LIBINT_SHGSHELL_ORDERING ${LIBINT_SHGSHELL_ORDERING_STANDARD})
endif()
if (LIBINT2_SHGAUSS_ORDERING STREQUAL "gaussian")
    set(LIBINT_SHGSHELL_ORDERING ${LIBINT_SHGSHELL_ORDERING_GAUSSIAN})
endif()

check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
if (HAVE_POSIX_MEMALIGN)
    set(LIBINT2_ALIGN_SIZE "0" CACHE STRING "(EXPERT) if posix_memalign is available, this will specify alignment of Libint data, in units of sizeof(LIBINT2_REALTYPE). Default is to use built-in heuristics")
    mark_as_advanced(LIBINT2_ALIGN_SIZE)
endif(HAVE_POSIX_MEMALIGN)

# Features =============================================================================================================

if (ENABLE_MPFR)
    find_package(MPFR)
    if (TARGET MPFR::GMPXX)
        set(LIBINT_HAS_MPFR 1)
    endif()
endif()

# look for Eigen
# prefer CMake-configured-and-installed instance
# re:NO_CMAKE_PACKAGE_REGISTRY: eigen3 registers its *build* tree with the user package registry ...
#                               to avoid issues with wiped build directory look for installed eigen
find_package(Eigen3 NO_MODULE QUIET NO_CMAKE_PACKAGE_REGISTRY)
if (TARGET Eigen3::Eigen)
  add_library(Eigen INTERFACE)
  foreach(prop INTERFACE_INCLUDE_DIRECTORIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS INTERFACE_LINK_LIBRARIES INTERFACE_POSITION_INDEPENDENT_CODE)
    get_property(EIGEN3_${prop} TARGET Eigen3::Eigen PROPERTY ${prop})
    set_property(TARGET Eigen PROPERTY
        ${prop} ${EIGEN3_${prop}})
  endforeach()
else (TARGET Eigen3::Eigen)
  # otherwise use FindEigen3.cmake module possibly installed somewhere in the path
  # but make sure EIGEN3_INCLUDE_DIR exists!
  find_package(Eigen3)
  if (EIGEN3_FOUND)
    if (NOT EXISTS "${EIGEN3_INCLUDE_DIR}")
      message(WARNING "Eigen3 is \"found\", but the reported EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR} does not exist; likely corrupt Eigen3 build registered in user or system package registry; specify EIGEN3_INCLUDE_DIR manually or (better) configure (with CMake) and install Eigen3 package")
    else(NOT EXISTS "${EIGEN3_INCLUDE_DIR}")
      add_library(Eigen INTERFACE)
      set_property(TARGET Eigen PROPERTY
          INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR})
    endif(NOT EXISTS "${EIGEN3_INCLUDE_DIR}")
  endif (EIGEN3_FOUND)
endif (TARGET Eigen3::Eigen)
if (TARGET Eigen)
  get_property(EIGEN3_INCLUDE_DIR TARGET Eigen PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Found Eigen3: EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}")
  install(TARGETS Eigen  EXPORT libint2 COMPONENT Eigen)
  set(LIBINT_HAS_CXX_API 1)
else(TARGET Eigen)
  set(LIBINT_HAS_CXX_API 0)
endif(TARGET Eigen)
if (REQUIRE_CXX_API AND NOT ${LIBINT_HAS_CXX_API})
  message(FATAL_ERROR "C++ API cannot be built without Eigen3; configure (via CMake) and install Eigen3 and add the install prefix to CMAKE_PREFIX_PATH, or add -DREQUIRE_CXX_API=OFF to the CMake command line if the C++ API is not required")
endif()

find_package(Boost 1.29)
if (Boost_FOUND)
    cmake_push_check_state()

    list(APPEND CMAKE_REQUIRED_FLAGS "-std=c++11")
    list(APPEND CMAKE_REQUIRED_INCLUDES ${Boost_INCLUDE_DIR})

    check_cxx_source_compiles("
        #include <boost/preprocessor.hpp>
        #if not BOOST_PP_VARIADICS  // no variadic macros? your compiler is out of date! (should not be possible since variadic macros are part of C++11)
        #  error \"your compiler does not provide variadic macros (but does support C++11), something is seriously broken, please create an issue at https://github.com/evaleev/libint/issues\"
        #endif
        int main() { return 0;}
        "
            LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS)

    cmake_pop_check_state()
else(Boost_FOUND)
    set(LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS 0)
endif(Boost_FOUND)

# if usable Boost.Preprocessor not found install the bundled version
if (NOT LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include/libint2)
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_SOURCE_DIR}/external/boost.tar.gz
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/include/libint2
            RESULT_VARIABLE UNPACK_BOOST_RESULT
            OUTPUT_VARIABLE UNPACK_BOOST_OUTPUT
            ERROR_VARIABLE UNPACK_BOOST_OUTPUT
    )
    if (NOT UNPACK_BOOST_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to unpack the bundled Boost! The tar command output:\n${UNPACK_BOOST_OUTPUT}")
    endif()
endif()

# Python is optional, unless ...
if (ENABLE_FORTRAN)  # ... need fortran
    find_package(PythonInterp REQUIRED)
else(ENABLE_FORTRAN)
    find_package(PythonInterp)
endif(ENABLE_FORTRAN)

# Set install paths ====================================================================================================

set(LIBINT2_INSTALL_BINDIR "bin"
        CACHE PATH "LIBINT2 binary install directory")
set(LIBINT2_INSTALL_INCLUDEDIR "include"
        CACHE PATH "LIBINT2 INCLUDE install directory")
set(LIBINT2_INSTALL_LIBDIR "lib"
        CACHE PATH "LIBINT2 LIB install directory")
set(LIBINT2_INSTALL_DATADIR "share/libint/${LIBINT_VERSION}"
        CACHE PATH "LIBINT2 DATA install directory")
set(LIBINT2_INSTALL_DOCDIR "share/libint/${LIBINT_VERSION}/doc"
        CACHE PATH "LIBINT2 DOC install directory")
set(LIBINT2_INSTALL_CMAKEDIR "lib/cmake/libint2"
        CACHE PATH "LIBINT2 CMAKE install directory")

# Libint library =======================================================================================================

include(srclist.cmake)
set(LIB_CXX_SRC )
foreach(FN IN LISTS LIBINT2_LIBRARY_CXX_SRC)
    list(APPEND LIB_CXX_SRC "src/${FN}")
endforeach()
# Create object files to use for static and shared libraries
add_library(int2_obj OBJECT ${LIB_CXX_SRC})
target_include_directories(int2_obj PRIVATE include ${PROJECT_BINARY_DIR}/include)
# Compile static library with position independent code

target_compile_definitions(int2_obj PRIVATE __COMPILING_LIBINT2)
target_compile_features(int2_obj PUBLIC "cxx_std_11")
if (TARGET MPFR::GMPXX)
  target_link_libraries(int2_obj PUBLIC MPFR::GMPXX)
endif()
if (NOT CMAKE_CXX_EXTENSIONS)
  set_target_properties(int2_obj PROPERTIES CXX_EXTENSIONS OFF)
endif(NOT CMAKE_CXX_EXTENSIONS)

if(BUILD_SHARED_LIBS OR LIBINT2_BUILD_SHARED_AND_STATIC_LIBS)
  set_target_properties(int2_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# shared and static libraries built from the same object files
if (LIBINT2_BUILD_SHARED_AND_STATIC_LIBS)
  add_library(int2 SHARED $<TARGET_OBJECTS:int2_obj>)
  add_library(int2_static STATIC $<TARGET_OBJECTS:int2_obj>)
  target_include_directories(int2_static INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/libint2>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${LIBINT2_INSTALL_INCLUDEDIR}>
        )
  target_compile_features(int2_static INTERFACE "cxx_std_11")
  # Add libraries to the list of installed components
  install(TARGETS int2_static EXPORT libint2
        COMPONENT int2
        LIBRARY DESTINATION "${LIBINT2_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${LIBINT2_INSTALL_LIBDIR}"
        # includes are installed by the include/CMakeLists.txt.include.export
        # INCLUDES DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}"
       )
else()
  add_library(int2 $<TARGET_OBJECTS:int2_obj>)
endif()
  # In any case, there will exist a target named "int2": diminish code length by
  # add just once the needed includes, flags and install directives.
  target_include_directories(int2 INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/libint2>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${LIBINT2_INSTALL_INCLUDEDIR}>
        )
  target_compile_features(int2 INTERFACE "cxx_std_11")
  # Add libraries to the list of installed components
  install(TARGETS int2 EXPORT libint2
        COMPONENT int2
        LIBRARY DESTINATION "${LIBINT2_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${LIBINT2_INSTALL_LIBDIR}"
        # includes are installed by the include/CMakeLists.txt.include.export
        # INCLUDES DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}"
       )

add_subdirectory(include)

if (LIBINT2_BUILD_SHARED_AND_STATIC_LIBS)
  set(int2_library int2_static)
else()
  set(int2_library int2)
endif()

# install basis set library
install(DIRECTORY ${CMAKE_SOURCE_DIR}/lib/basis
        COMPONENT int2
        DESTINATION "${LIBINT2_INSTALL_DATADIR}"
        )

# LibintCXX library ====================================================================================================

if (LIBINT_HAS_CXX_API)
  add_library(cxx INTERFACE)
  target_compile_features(cxx INTERFACE "cxx_std_11")
  target_link_libraries(cxx INTERFACE ${int2_library} Eigen)
  if (LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS)
      target_include_directories(cxx INTERFACE ${Boost_INCLUDE_DIR})
  endif(LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS)
  get_filename_component(DATADIR_ABSOLUTE "${CMAKE_INSTALL_PREFIX}/${LIBINT2_INSTALL_DATADIR}" ABSOLUTE)
  target_compile_definitions(cxx INTERFACE
          $<BUILD_INTERFACE:SRCDATADIR="${CMAKE_SOURCE_DIR}/lib/basis">
          $<INSTALL_INTERFACE:DATADIR="${DATADIR_ABSOLUTE}">)
  # Add library to the list of installed components
  install(TARGETS cxx EXPORT libint2
          COMPONENT cxx
          LIBRARY DESTINATION "${LIBINT2_INSTALL_LIBDIR}"
          ARCHIVE DESTINATION "${LIBINT2_INSTALL_LIBDIR}"
          # includes are installed by the include/CMakeLists.txt.include.export
          # INCLUDES DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}"
          )
endif(LIBINT_HAS_CXX_API)

# Tests ================================================================================================================

enable_testing(true)
add_custom_target(check USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} -V)

add_executable(eritest EXCLUDE_FROM_ALL tests/eri/test.cc)
target_link_libraries(eritest ${int2_library})
target_compile_definitions(eritest PRIVATE __COMPILING_LIBINT2)
target_include_directories(eritest PRIVATE ${CMAKE_SOURCE_DIR}/tests/eri)

add_test(eritest/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target eritest)
set_tests_properties(eritest/build PROPERTIES FIXTURES_SETUP ERITEST_EXEC)
add_test(NAME eritest/run0
        COMMAND $<TARGET_FILE:eritest> 0 2)
set_tests_properties(eritest/run0
        PROPERTIES FIXTURES_REQUIRED ERITEST_EXEC)
add_test(NAME eritest/run1
        COMMAND $<TARGET_FILE:eritest> 1 1)
set_tests_properties(eritest/run1
        PROPERTIES FIXTURES_REQUIRED ERITEST_EXEC)
add_test(NAME eritest/run2
        COMMAND $<TARGET_FILE:eritest> 2 1)
set_tests_properties(eritest/run2
        PROPERTIES FIXTURES_REQUIRED ERITEST_EXEC)

if (LIBINT_HAS_CXX_API)
    set(utests_src
            tests/unit/c-api.c
            tests/unit/c-api-util.cc
            tests/unit/catch.hpp
            tests/unit/fixture.h
            tests/unit/test-1body.cc
            tests/unit/test-2body.cc
            tests/unit/test-core-ints.cc
            tests/unit/test-c-api.cc
            tests/unit/test-core.cc
            tests/unit/test-permute.cc
            tests/unit/test-shell-order.cc
            tests/unit/test-util.cc
            )
    add_executable(unit_tests EXCLUDE_FROM_ALL tests/unit/test.cc ${utests_src})
    target_compile_definitions(unit_tests PRIVATE __COMPILING_LIBINT2)
    target_link_libraries(unit_tests cxx)
    add_test(unit/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target unit_tests)
    set_tests_properties(unit/build PROPERTIES FIXTURES_SETUP UNIT_TESTS_EXEC)
    add_test(NAME unit/run
            COMMAND $<TARGET_FILE:unit_tests>)
    set_tests_properties(unit/run
            PROPERTIES FIXTURES_REQUIRED UNIT_TESTS_EXEC)

    add_executable(hartree-fock EXCLUDE_FROM_ALL tests/hartree-fock/hartree-fock.cc)
    target_compile_definitions(hartree-fock PRIVATE __COMPILING_LIBINT2)
    target_link_libraries(hartree-fock cxx)
    add_test(hftest/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hartree-fock)
    set_tests_properties(hftest/build PROPERTIES FIXTURES_SETUP HFTEST_EXEC)
    if (PYTHONINTERP_FOUND)
        add_test(NAME hftest/run
                COMMAND COMMAND ${CMAKE_COMMAND}
                -DtestName=hartree-fock
                -DtestArgs=${CMAKE_SOURCE_DIR}/tests/hartree-fock/h2o.xyz
                -DsrcDir=${CMAKE_SOURCE_DIR}
                -DpythonExec=${PYTHON_EXECUTABLE}
                -P ${CMAKE_SOURCE_DIR}/cmake/hftest.cmake)
    else()
        add_test(NAME hftest/run
                COMMAND $<TARGET_FILE:hartree-fock> ${CMAKE_SOURCE_DIR}/tests/hartree-fock/h2o.xyz)
    endif()
    set_tests_properties(hftest/run
            PROPERTIES FIXTURES_REQUIRED HFTEST_EXEC)

    find_package(Threads)  # for some reason clang does not link in threading support even though we are using C++ threads
    add_executable(hartree-fock++ EXCLUDE_FROM_ALL tests/hartree-fock/hartree-fock++.cc)
    target_compile_definitions(hartree-fock++ PRIVATE __COMPILING_LIBINT2)
    target_link_libraries(hartree-fock++ cxx ${CMAKE_THREAD_LIBS_INIT})
    add_test(hf++test/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hartree-fock++)
    set_tests_properties(hf++test/build PROPERTIES FIXTURES_SETUP HFXXTEST_EXEC)
    if (PYTHONINTERP_FOUND)
        add_test(NAME hf++test/run
                COMMAND COMMAND ${CMAKE_COMMAND}
                -DtestName=hartree-fock++
                -DtestArgs=${CMAKE_SOURCE_DIR}/tests/hartree-fock/h2o_rotated.xyz
                -DsrcDir=${CMAKE_SOURCE_DIR}
                -DpythonExec=${PYTHON_EXECUTABLE}
                -P ${CMAKE_SOURCE_DIR}/cmake/hftest.cmake)
    else()
        add_test(NAME hf++test/run
                COMMAND $<TARGET_FILE:hartree-fock++> ${CMAKE_SOURCE_DIR}/tests/hartree-fock/h2o_rotated.xyz)
    endif()
    set_tests_properties(hf++test/run
            PROPERTIES FIXTURES_REQUIRED HFXXTEST_EXEC)
endif (LIBINT_HAS_CXX_API)

# Fortran bindings =====================================================================================================

if (ENABLE_FORTRAN)
    # preprocess libint2.h ... this is a guess for non-unix
    if (UNIX)

        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)

        # preprocessed libint.h
        add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/fortran/libint2.h.i
                COMMAND ${CMAKE_C_COMPILER} -E -D__COMPILING_LIBINT2 -I${PROJECT_SOURCE_DIR}/include -I${PROJECT_BINARY_DIR}/include
                -I${PROJECT_BINARY_DIR}/include/libint2 ${CMAKE_SOURCE_DIR}/include/libint2.h -o ${CMAKE_BINARY_DIR}/fortran/libint2.h.i
                DEPENDS ${CMAKE_SOURCE_DIR}/include/libint2.h
                COMMENT "Generating libint2.h.i"
                )
    else()
        message(FATAL_ERROR "Cannot run preprocessor on non-Unix systems, disable Fortran to proceed")
    endif()

    # translated Libint_t
    add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/fortran/libint2_types_f.h
            COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/fortran/c_to_f.py ${CMAKE_BINARY_DIR}/fortran/libint2.h.i ${CMAKE_BINARY_DIR}/fortran/libint2_types_f.h Libint_t
            DEPENDS ${CMAKE_BINARY_DIR}/fortran/libint2.h.i
            COMMENT "Generating libint2_types_f.h"
            )

    # extracted defines from libint2_types.h
    add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/fortran/fortran_incldefs.h
            COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/fortran/make_defs.py ${CMAKE_SOURCE_DIR}/include/libint2_types.h ${CMAKE_BINARY_DIR}/fortran/fortran_incldefs.h
            DEPENDS ${CMAKE_SOURCE_DIR}/include/libint2_types.h
            COMMENT "Generating fortran_incldefs.h"
            )

    # build module
    add_library(libint_f OBJECT fortran/libint_f.F90)
    set_source_files_properties(fortran/libint_f.F90 PROPERTIES OBJECT_DEPENDS "${CMAKE_BINARY_DIR}/fortran/libint2_types_f.h;${CMAKE_BINARY_DIR}/fortran/fortran_incldefs.h")
    target_include_directories(libint_f PUBLIC
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/fortran>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${LIBINT2_INSTALL_INCLUDEDIR}>)
    target_compile_definitions(libint_f PRIVATE __COMPILING_LIBINT2)
    add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/libint_f.mod
                       COMMAND ${CMAKE_COMMAND} -E "message(STATUS \"Building Fortran module libint_f.mod\")"
                       DEPENDS libint_f
                       VERBATIM USES_TERMINAL)

    # tests
    add_executable(fortran_example EXCLUDE_FROM_ALL fortran/fortran_example.F90 $<TARGET_OBJECTS:libint_f>)
    target_compile_definitions(fortran_example PRIVATE __COMPILING_LIBINT2)
    target_link_libraries(fortran_example ${int2_library})

    add_test(fortran_example/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target fortran_example)
    set_tests_properties(fortran_example/build PROPERTIES FIXTURES_SETUP FORTRAN_EXAMPLE_EXEC)
    add_test(NAME fortran_example/run
            COMMAND $<TARGET_FILE:fortran_example>)
    set_tests_properties(fortran_example/run
            PROPERTIES FIXTURES_REQUIRED FORTRAN_EXAMPLE_EXEC)

    if (LIBINT_HAS_CXX_API)
        add_executable(fortran_test EXCLUDE_FROM_ALL fortran/test.cc fortran/test-eri.cc $<TARGET_OBJECTS:libint_f>)
        target_compile_definitions(fortran_test PRIVATE __COMPILING_LIBINT2)
        target_link_libraries(fortran_test cxx)
        add_test(fortran_test/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target fortran_test)
        set_tests_properties(fortran_test/build PROPERTIES FIXTURES_SETUP FORTRAN_TEST_EXEC)
        add_test(NAME fortran_test/run
                COMMAND $<TARGET_FILE:fortran_test>)
        set_tests_properties(fortran_test/run
                PROPERTIES FIXTURES_REQUIRED FORTRAN_TEST_EXEC)
    endif (LIBINT_HAS_CXX_API)

    # install targets
    install(FILES ${CMAKE_BINARY_DIR}/libint_f.mod
            COMPONENT int2
            DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}"
            )

endif(ENABLE_FORTRAN)

# Configure files ======================================================================================================

configure_file(
        ${PROJECT_SOURCE_DIR}/include/libint2/config.h.cmake.in
        ${PROJECT_BINARY_DIR}/include/libint2/config.h
        @ONLY
)

configure_file(
        include/libint2/basis.h.in
        ${PROJECT_BINARY_DIR}/include/libint2/basis.h
        @ONLY
)

configure_file(
        libint2.pc.in
        ${PROJECT_BINARY_DIR}/libint2.pc
        @ONLY
)

# install cmake-processed files
install(FILES ${PROJECT_BINARY_DIR}/libint2.pc
        DESTINATION lib/pkgconfig)
install(FILES ${PROJECT_BINARY_DIR}/include/libint2/config.h
        DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}/libint2")
install(FILES ${PROJECT_BINARY_DIR}/include/libint2/basis.h
        DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}/libint2")

# Create the version file
write_basic_package_version_file(libint2-config-version.cmake
        VERSION ${LIBINT_VERSION} COMPATIBILITY AnyNewerVersion)

## Create the configure file
configure_package_config_file(cmake/libint2-config.cmake.in
        "${PROJECT_BINARY_DIR}/libint2-config.cmake"
        INSTALL_DESTINATION "${LIBINT2_INSTALL_CMAKEDIR}"
        PATH_VARS CMAKE_INSTALL_PREFIX LIBINT2_INSTALL_BINDIR
        LIBINT2_INSTALL_INCLUDEDIR LIBINT2_INSTALL_LIBDIR
        LIBINT2_INSTALL_DOCDIR LIBINT2_INSTALL_CMAKEDIR)

## Install config, version, and target files
install(EXPORT libint2
        FILE "libint2-targets.cmake"
        DESTINATION "${LIBINT2_INSTALL_CMAKEDIR}"
        NAMESPACE Libint2::
        COMPONENT libint2-config)
install(FILES
        "${PROJECT_BINARY_DIR}/libint2-config.cmake"
        "${PROJECT_BINARY_DIR}/libint2-config-version.cmake"
        DESTINATION "${LIBINT2_INSTALL_CMAKEDIR}"
        COMPONENT libint2-config)
add_custom_target(install-config
        COMMAND ${CMAKE_COMMAND} -DCOMPONENT=libint2-config -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
        COMMENT "Installing Libint2 config components")

feature_summary(WHAT ALL
        DESCRIPTION "=== Libint2 Package/Feature Info ===")

###############################################################################
# appendix: misc details
###############################################################################
SET(CMAKE_COLOR_MAKEFILE ON)
