list(APPEND CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR}/cmake/modules)

include(GNUInstallDirs)
include(int_computed)  # for macros.tex

# <<<  Generate Library Source  >>>

file(MAKE_DIRECTORY "${EXPORT_STAGE_DIR}/src")
add_custom_command(OUTPUT "${EXPORT_STAGE_DIR}/src/libint2_params.h"
        COMMAND "${PROJECT_BINARY_DIR}/src/bin/libint/build_libint"
        WORKING_DIRECTORY "${EXPORT_STAGE_DIR}/src"
        DEPENDS "${PROJECT_BINARY_DIR}/src/bin/libint/build_libint"
        COMMENT "Generating Libint2 library source")
add_custom_target(libint-library-generate DEPENDS "${EXPORT_STAGE_DIR}/src/libint2_params.h")

# <<<  Add Metadata To The Library Source  >>>

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
        -P ${CMAKE_CURRENT_SOURCE_DIR}/populate.cmake
        COMMENT "Populating Libint2 library")
add_custom_target(libint-library-populate DEPENDS ${EXPORT_STAGE_DIR}/CMakeLists.txt)

# <<<  Export The Library Source  >>>

add_custom_command(OUTPUT "${EXPORT_STAGE_DIR}.tgz"
        COMMAND ${CMAKE_COMMAND} -E tar "cfz" "${EXPORT_STAGE_DIR}.tgz" "${EXPORT_STAGE_DIR}"
        WORKING_DIRECTORY "${EXPORT_STAGE_DIR}/.."
        DEPENDS libint-library-generate libint-library-populate
        COMMENT "Exporting tarball of Libint2 library source")
if (LIBINT_BUILD_LIBRARY_AS_SUBPROJECT)
    # rerun cmake if building as subproject to unpack and consume the library as subproject
    add_custom_target(libint-library-export DEPENDS "${EXPORT_STAGE_DIR}.tgz"
            COMMAND ${CMAKE_COMMAND} -S "${CMAKE_SOURCE_DIR}" -B "${CMAKE_BINARY_DIR}")
else()
    add_custom_target(libint-library-export DEPENDS "${EXPORT_STAGE_DIR}.tgz")
endif()
add_custom_target(export DEPENDS libint-library-export)


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

else()

    include(ExternalProject)
    include(GNUInstallDirs)

    list(APPEND library_CMAKE_ARGS
            -G${CMAKE_GENERATOR}
            -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
            -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
            -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
            -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
            -DLIBINT2_INSTALL_CMAKEDIR=${LIBINT2_INSTALL_CMAKEDIR}
            -DLIBINT2_INSTALL_BASISDIR=${LIBINT2_INSTALL_BASISDIR}
            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            -DLIBINT2_BUILD_SHARED_AND_STATIC_LIBS=${LIBINT2_BUILD_SHARED_AND_STATIC_LIBS}
            -DREQUIRE_CXX_API=${REQUIRE_CXX_API}
            -DREQUIRE_CXX_API_COMPILED=${REQUIRE_CXX_API_COMPILED}
            -DENABLE_FORTRAN=${ENABLE_FORTRAN}
            -DENABLE_PYTHON=${ENABLE_PYTHON}
            -DPREFIX_PYTHON_INSTALL=${PREFIX_PYTHON_INSTALL}
            -DLIBINT2_SHGAUSS_ORDERING=${LIBINT2_SHGAUSS_ORDERING}
            -DBUILD_TESTING=${BUILD_TESTING}
            -DENABLE_MPFR=${ENABLE_MPFR}
            -DLIBINT2_REALTYPE=${LIBINT2_REALTYPE}
            -DENABLE_XHOST=${ENABLE_XHOST}
            )
    if (CMAKE_C_COMPILER)
        list(APPEND library_CMAKE_ARGS
                "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
                )
    endif()
    if (ENABLE_FORTRAN AND CMAKE_Fortran_COMPILER)
        list(APPEND library_CMAKE_ARGS
                "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
                )
    endif()
    if (Eigen3_DIR)
        list(APPEND library_CMAKE_ARGS
                "-DEigen3_DIR=${Eigen3_DIR}"
                )
    endif()
    if (Eigen3_ROOT)
        list(APPEND library_CMAKE_ARGS
                "-DEigen3_ROOT=${Eigen3_ROOT}"
                )
    endif()
    if (Boost_DIR)
        list(APPEND library_CMAKE_ARGS
                "-DBoost_DIR=${Boost_DIR}"
                )
    endif()
    if (BOOST_ROOT)
        list(APPEND library_CMAKE_ARGS
                "-DBOOST_ROOT=${BOOST_ROOT}"
                )
    endif()
#    if (CMAKE_DISABLE_FIND_PACKAGE_Boost)
#        list(APPEND library_CMAKE_ARGS
#                "-DCMAKE_DISABLE_FIND_PACKAGE_Boost=${CMAKE_DISABLE_FIND_PACKAGE_Boost}"
#    endif()
    if (LIBINT_LOCAL_Eigen3_INSTALL)
        list(APPEND library_CMAKE_ARGS
                "-DLIBINT_LOCAL_Eigen3_INSTALL=${LIBINT_LOCAL_Eigen3_INSTALL}"
                )
    endif()
    message(STATUS "invoking ExternalProject_add to build+test Libint library with library_CMAKE_ARGS=${library_CMAKE_ARGS}")

    ExternalProject_Add(library
            DEPENDS libint-library-export
            BUILD_ALWAYS 1
            # download generalized to generating source by unpacking source tarball
            URL "${EXPORT_STAGE_DIR}.tgz"
            # DOWNLOAD_DIR <SOURCE_DIR>  # this is tidier but fails on Windows
            LOG_DOWNLOAD 1
            CMAKE_ARGS "${library_CMAKE_ARGS}"
            CMAKE_CACHE_ARGS
              -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
              -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
              -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH}
            TEST_COMMAND ""
            USES_TERMINAL_CONFIGURE TRUE
            USES_TERMINAL_BUILD TRUE
            USES_TERMINAL_INSTALL TRUE
            STEP_TARGETS configure
            BUILD_BYPRODUCTS ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR} ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR} ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}
            )
    ExternalProject_Get_property(library BINARY_DIR)

    # instead of running tests by default, replicate check-libint2 target here and use it to invoke library's check-libint2 target
    if (BUILD_TESTING)
      include(AddCustomTargetSubproject)
      add_custom_target_subproject(libint2 check USES_TERMINAL
              COMMAND ${CMAKE_CTEST_COMMAND} -V -R "libint2/"
              WORKING_DIRECTORY ${BINARY_DIR}
              DEPENDS library-configure)
    endif()

    install(
      DIRECTORY ${STAGED_INSTALL_PREFIX}/
      DESTINATION ${CMAKE_INSTALL_PREFIX}
      USE_SOURCE_PERMISSIONS OPTIONAL
      )

endif()
