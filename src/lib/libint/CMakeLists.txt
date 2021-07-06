cmake_minimum_required(VERSION 3.8)
cmake_policy(SET CMP0074 NEW)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR}/cmake/modules)

include(GNUInstallDirs)
include(int_computed)  # for macros.tex

# <<<  Generate Library  >>>

file(MAKE_DIRECTORY "${EXPORT_STAGE_DIR}/src")
add_custom_command(OUTPUT "${EXPORT_STAGE_DIR}/src/libint2_params.h"
        COMMAND "${PROJECT_BINARY_DIR}/src/bin/libint/build_libint"
        WORKING_DIRECTORY "${EXPORT_STAGE_DIR}/src"
        DEPENDS "${PROJECT_BINARY_DIR}/src/bin/libint/build_libint"
        COMMENT "Generating Libint2 library")
add_custom_target(libint-library-generate DEPENDS "${EXPORT_STAGE_DIR}/src/libint2_params.h")

# <<<  Add Metadata To The Library >>>

configure_file("${PROJECT_SOURCE_DIR}/doc/progman/macros.tex.in" "${EXPORT_STAGE_DIR}/doc/macros.tex" @ONLY)
set(DATADIR_ABSOLUTE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION})
configure_file("${PROJECT_SOURCE_DIR}/include/libint2/basis.h.in" "${EXPORT_STAGE_DIR}/include/libint2/basis.h" @ONLY)
configure_file("features.in" "${EXPORT_STAGE_DIR}/features" @ONLY)

add_custom_command(OUTPUT ${EXPORT_STAGE_DIR}/CMakeLists.txt
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/populate.cmake"
                libint-library-generate
        COMMAND ${CMAKE_COMMAND}
        "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
        "-DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}"
        "-DEXPORT_STAGE_DIR=${EXPORT_STAGE_DIR}"
        "-DLIBINT_VERSION=${LIBINT_VERSION}"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/populate.cmake
        COMMENT "Populating Libint2 library")
add_custom_target(libint-library-populate DEPENDS ${EXPORT_STAGE_DIR}/CMakeLists.txt)

# <<<  Export The Library >>>

add_custom_command(OUTPUT "${EXPORT_STAGE_DIR}.tgz"
        COMMAND ${CMAKE_COMMAND} -E tar "cfvz" "${EXPORT_STAGE_DIR}.tgz" "${EXPORT_STAGE_DIR}"
        WORKING_DIRECTORY "${EXPORT_STAGE_DIR}/.."
        DEPENDS libint-library-generate libint-library-populate
        COMMENT "Exporting Libint2 library")
if (LIBINT_BUILD_LIBRARY_AS_SUBPROJECT)
    # rerun cmake if building as subproject to unpack and consume the library as subproject
    add_custom_target(libint-library-export DEPENDS "${EXPORT_STAGE_DIR}.tgz"
            COMMAND ${CMAKE_COMMAND} -S "${CMAKE_SOURCE_DIR}" -B "${CMAKE_BINARY_DIR}")
else (LIBINT_BUILD_LIBRARY_AS_SUBPROJECT)
    add_custom_target(libint-library-export DEPENDS "${EXPORT_STAGE_DIR}.tgz")
endif (LIBINT_BUILD_LIBRARY_AS_SUBPROJECT)


    # <<<  Build Library  >>>

# make sure CMAKE_C_COMPILER is defined
include(CheckLanguage)
check_language(C)
if(CMAKE_C_COMPILER)
    enable_language(C)
else()
    message(FATAL_ERROR "No C language support; invoke CMake with -DCMAKE_C_COMPILER=/path/to/C/compiler to resolve")
endif()

if (LIBINT_BUILD_LIBRARY_AS_SUBPROJECT)

    if (EXISTS "${EXPORT_STAGE_DIR}.tgz")
        message("Unpacking ${EXPORT_STAGE_DIR}.tgz and setting up as a subproject")
        include(FetchContent)
        FetchContent_Declare(
            libint-library
            URL "${EXPORT_STAGE_DIR}.tgz"
        )
        FetchContent_MakeAvailable(libint-library)
    else()
        message(WARNING "LIBINT_BUILD_LIBRARY_AS_SUBPROJECT=ON but library has not been exported yet; build target 'libint-library-export' first, this will make library build/test targets ('check', 'install') available")
    endif()

else(LIBINT_BUILD_LIBRARY_AS_SUBPROJECT)

    include(ExternalProject)
    include(GNUInstallDirs)

    ExternalProject_Add(library
            DEPENDS libint-library-export
            BUILD_ALWAYS 1
            # download generalized to generating source by unpacking source tarball
            URL "${EXPORT_STAGE_DIR}.tgz"
            # DOWNLOAD_DIR <SOURCE_DIR>  # this is tidier but fails on Windows
            LOG_DOWNLOAD 1
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}  # for a C API test
            -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
            -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
            -DBUILD_SHARED=${BUILD_SHARED}
            -DBUILD_STATIC=${BUILD_STATIC}
            -DENABLE_XHOST=${ENABLE_XHOST}
            -DBUILD_FPIC=${BUILD_FPIC}
            -DBOOST_ROOT=${BOOST_ROOT}
            -DEigen3_ROOT=${Eigen3_ROOT}
            -DMPFR_ROOT=${MPFR_ROOT}
            -DLIBINT2_SHGAUSS_ORDERING=${LIBINT2_SHGAUSS_ORDERING}
            -DBUILD_TESTING=${BUILD_TESTING}
            CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
            TEST_COMMAND ${CMAKE_CTEST_COMMAND} --extra-verbose --output-on-failure)

    install(DIRECTORY ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
            # BINDIR excluded b/c `build_libint` not really for distribution
            ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
            ${STAGED_INSTALL_PREFIX}/share
            DESTINATION ${CMAKE_INSTALL_PREFIX}
            USE_SOURCE_PERMISSIONS)

endif(LIBINT_BUILD_LIBRARY_AS_SUBPROJECT)
