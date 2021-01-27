include(ExternalProject)
include(GNUInstallDirs)

add_subdirectory(bin/libint)
add_subdirectory(lib/libint)

# <<<  Build Library  >>>

ExternalProject_Add(library
    DEPENDS libint-library-export
    BUILD_ALWAYS 1
    # download generalized to generating source by unpacking source tarball
    URL "${EXPORT_STAGE_DIR}.tgz"
    DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/library-build"
    # DOWNLOAD_DIR <SOURCE_DIR>  # this is tidier but fails on Windows
    LOG_DOWNLOAD 1
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}"
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
