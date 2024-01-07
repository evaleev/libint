list(APPEND CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR}/cmake/modules)

include(GNUInstallDirs)
include(int_computed)  # for macros.tex and features

# <<<  Generate Library Source  >>>

file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/generated")
add_custom_command(
  OUTPUT
    "${PROJECT_BINARY_DIR}/generated/libint2_params.h"
    "${PROJECT_BINARY_DIR}/generated/HRRPart0bra0ket0pp.h"
    "${PROJECT_BINARY_DIR}/generated/HRRPart0bra0ket0pp.cc"
  COMMAND $<TARGET_FILE:build_libint>
  DEPENDS
    build_libint
  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/generated"
  COMMENT "Generating Libint2 library source"
  )
add_custom_target(libint-library-generate DEPENDS "${PROJECT_BINARY_DIR}/generated/libint2_params.h")


# <<<  Add Metadata To The Library Source  >>>

# the library source to be tarballed is mainly assembled through populate.cmake.
#   These following configure_file's need a few extra substitution variables that
#   aren't passed to populate.

configure_file(
  "features.cmake.in"
  "${EXPORT_STAGE_DIR}/features"
  @ONLY)

configure_file(
  "${PROJECT_SOURCE_DIR}/doc/progman/macros.tex.in"
  "${EXPORT_STAGE_DIR}/doc/macros.tex"
  @ONLY)

# position basis.h for local buildtree.
# * gets re-configured (CMakeLists.txt.export) for library export with its own
#   CMAKE_INSTALL_PREFIX
# * gets re-re-positioned (end of this file) for case of continuous generator +
#   library build (through EP) with outer CMAKE_INSTALL_PREFIX
set(DATADIR_ABSOLUTE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION})
configure_file(
  "${PROJECT_SOURCE_DIR}/include/libint2/basis.h.in"  # TODO basis.h.cmake.in after libtool retires
  "include/libint2/basis.h"
  @ONLY)

add_custom_command(
  OUTPUT
    ${EXPORT_STAGE_DIR}/CMakeLists.txt
    ${EXPORT_STAGE_DIR}/include/libint2_params.h
    ${EXPORT_STAGE_DIR}/include/HRRPart0bra0ket0pp.h
    ${EXPORT_STAGE_DIR}/src/HRRPart0bra0ket0pp.cc
  COMMAND ${CMAKE_COMMAND}
          "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
          "-DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}"
          "-DEXPORT_STAGE_DIR=${EXPORT_STAGE_DIR}"
          -P ${CMAKE_CURRENT_SOURCE_DIR}/populate.cmake
  DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/populate.cmake"
    "${PROJECT_BINARY_DIR}/generated/libint2_params.h"
    "${PROJECT_BINARY_DIR}/generated/HRRPart0bra0ket0pp.h"
    "${PROJECT_BINARY_DIR}/generated/HRRPart0bra0ket0pp.cc"
    libint-library-generate
  COMMENT "Populating Libint2 library"
  )
add_custom_target(libint-library-populate DEPENDS ${EXPORT_STAGE_DIR}/CMakeLists.txt)


# <<<  Export The Library Source  >>>

set(EXPORT_TARBALL "${EXPORT_STAGE_DIR}.${EXPORT_COMPRESSOR_EXT}")
add_custom_command(
  OUTPUT
    "${EXPORT_TARBALL}"
  COMMAND ${CMAKE_COMMAND} -E tar ${EXPORT_COMPRESSOR_CMD} "${EXPORT_TARBALL}" "${EXPORT_STAGE_DIR}"
  DEPENDS
    libint-library-generate
    libint-library-populate
    ${EXPORT_STAGE_DIR}/CMakeLists.txt
  WORKING_DIRECTORY "${EXPORT_STAGE_DIR}/.."
  COMMENT "Exporting tarball of Libint2 library source: ${EXPORT_TARBALL}"
  )
if (LIBINT_BUILD_LIBRARY_AS_SUBPROJECT)
    # rerun cmake if building as subproject to unpack and consume the library as subproject
    add_custom_target(libint-library-export DEPENDS "${EXPORT_TARBALL}"
      COMMAND ${CMAKE_COMMAND} -S "${CMAKE_SOURCE_DIR}" -B "${CMAKE_BINARY_DIR}")
else()
    add_custom_target(libint-library-export DEPENDS "${EXPORT_TARBALL}")
endif()
add_custom_target(export DEPENDS libint-library-export)


# <<<  Build Library  >>>

include(CheckLanguage)
check_language(C)
if(CMAKE_C_COMPILER)
    enable_language(C)
else()
    message(FATAL_ERROR "No C language support; invoke CMake with -D CMAKE_C_COMPILER=/path/to/C/compiler to resolve")
endif()

if (LIBINT_BUILD_LIBRARY_AS_SUBPROJECT)
    include(FetchContent)

    if (EXISTS "${EXPORT_TARBALL}")
        message("Unpacking ${EXPORT_TARBALL} and setting up as a subproject")
        FetchContent_Declare(
          libint-library
          URL "${EXPORT_TARBALL}"
          )
        FetchContent_MakeAvailable(libint-library)
    else()
        message(WARNING "LIBINT_BUILD_LIBRARY_AS_SUBPROJECT=ON but library has not been exported yet; build target 'libint-library-export' first, this will make library build/test targets ('check', 'install') available. Don't run 'check' with CTEST_PARALLEL_LEVEL as this will build the library multiple times.")
    endif()

else()
    include(ExternalProject)

    list(APPEND library_CMAKE_ARGS
        -G${CMAKE_GENERATOR}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBUILD_TESTING=${BUILD_TESTING}
        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
        -DLIBINT2_BUILD_SHARED_AND_STATIC_LIBS=${LIBINT2_BUILD_SHARED_AND_STATIC_LIBS}
        -DLIBINT2_REQUIRE_CXX_API=${LIBINT2_REQUIRE_CXX_API}
        -DLIBINT2_REQUIRE_CXX_API_COMPILED=${LIBINT2_REQUIRE_CXX_API_COMPILED}
        -DLIBINT2_ENABLE_FORTRAN=${LIBINT2_ENABLE_FORTRAN}
        -DLIBINT2_ALIGN_SIZE=${LIBINT_ALIGN_SIZE}
        -DLIBINT2_REALTYPE=${LIBINT2_REALTYPE}
        -DLIBINT2_ENABLE_MPFR=${LIBINT2_ENABLE_MPFR}
        -DLIBINT2_ENABLE_PYTHON=${LIBINT2_ENABLE_PYTHON}
        -DLIBINT2_PREFIX_PYTHON_INSTALL=${LIBINT2_PREFIX_PYTHON_INSTALL}
      )
    if (Multiprecision_ROOT)
        list(APPEND library_CMAKE_ARGS "-DMultiprecision_ROOT=${Multiprecision_ROOT}")
    endif()
    if (Boost_DIR)
        list(APPEND library_CMAKE_ARGS "-DBoost_DIR=${Boost_DIR}")
    endif()
    if (BOOST_ROOT)
        list(APPEND library_CMAKE_ARGS "-DBOOST_ROOT=${BOOST_ROOT}")
    endif()
    if (CMAKE_DISABLE_FIND_PACKAGE_Boost)
        list(APPEND library_CMAKE_ARGS "-DCMAKE_DISABLE_FIND_PACKAGE_Boost=${CMAKE_DISABLE_FIND_PACKAGE_Boost}")
    endif()
    if (Eigen3_DIR)
        list(APPEND library_CMAKE_ARGS "-DEigen3_DIR=${Eigen3_DIR}")
    endif()
    if (Eigen3_ROOT)
        list(APPEND library_CMAKE_ARGS "-DEigen3_ROOT=${Eigen3_ROOT}")
    endif()
    if (LIBINT2_LOCAL_Eigen3_INSTALL)
        list(APPEND library_CMAKE_ARGS "-DLIBINT2_LOCAL_Eigen3_INSTALL=${LIBINT2_LOCAL_Eigen3_INSTALL}")
    endif()
    if (Python_EXECUTABLE)
        list(APPEND library_CMAKE_ARGS "-DPython_EXECUTABLE=${Python_EXECUTABLE}")
    endif()
    if (LIBINT2_ENABLE_FORTRAN AND CMAKE_Fortran_COMPILER)
        # Note: not handling CMAKE_Fortran_COMPILER_ARG1
        list(APPEND library_CMAKE_ARGS "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
    endif()

    ExternalProject_Add(
      library
      DEPENDS
        libint-library-export
      BUILD_ALWAYS 1
      # download generalized to generating source by unpacking source tarball
      # * DOWNLOAD_DIR <SOURCE_DIR> is tidier but fails on Windows
      URL "${EXPORT_TARBALL}"
      LOG_DOWNLOAD 1
      CMAKE_ARGS "${library_CMAKE_ARGS}"
      CMAKE_CACHE_ARGS
        -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
        -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
        -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH}
      TEST_COMMAND ""
      USES_TERMINAL_CONFIGURE TRUE
      USES_TERMINAL_BUILD TRUE
      USES_TERMINAL_INSTALL TRUE
      STEP_TARGETS
        build
        install
      )
    ExternalProject_Get_property(library BINARY_DIR)

    # instead of running tests by default, replicate check-libint2 target here
    #   and use it to invoke library's check-libint2 target
    if (BUILD_TESTING)
        include(AddCustomTargetSubproject)
        add_custom_target_subproject(
          libint2 check  # target check-libint2
          USES_TERMINAL
          COMMAND ${CMAKE_CTEST_COMMAND} -V -R "libint2/"
          WORKING_DIRECTORY ${BINARY_DIR}
          DEPENDS
            "library-build"
          )

        if (LIBINT2_ENABLE_PYTHON)
            find_package(Python COMPONENTS Interpreter REQUIRED)
            add_custom_target_subproject(
              libint2 check-python
              COMMAND ${Python_EXECUTABLE} -m setup test
              WORKING_DIRECTORY ${BINARY_DIR}/python
              DEPENDS
                "library-install"
              )
        endif()
    endif()

    install(
      FILES
        "${CMAKE_CURRENT_BINARY_DIR}/include/libint2/basis.h"
      DESTINATION "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/libint2"
      )
endif()
