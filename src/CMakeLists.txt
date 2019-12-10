include(ExternalProject)
include(GNUInstallDirs)

ExternalProject_Add(generator
    BUILD_ALWAYS 1
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin/libint
    CMAKE_ARGS -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
               -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DLibint2_BINARY_DIR=${Libint2_BINARY_DIR}
               -DLibint2_ROOT=${Libint2_SOURCE_DIR}
               -DBOOST_ROOT=${BOOST_ROOT}
               -DMPFR_ROOT=${MPFR_ROOT}
    CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    TEST_COMMAND ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/test
    TEST_AFTER_INSTALL TRUE)

ExternalProject_Add(export
    DEPENDS generator
    BUILD_ALWAYS 1
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/libint
    DOWNLOAD_DIR ${EXPORT_INSTALL_PREFIX}/src
    # download generalized to generating source by running build_libint generator
    DOWNLOAD_COMMAND ${STAGED_INSTALL_PREFIX}/bin/build_libint
    LOG_DOWNLOAD 1
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXPORT_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}  # unused but since CXX declared will need to detect compiler
               -DLibint2_BINARY_DIR=${Libint2_BINARY_DIR}
               -DLibint2_ROOT=${Libint2_SOURCE_DIR}
    CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS})
ExternalProject_Add_Step(export tar
    DEPENDEES install
    COMMAND ${CMAKE_COMMAND} -E tar "cfvz" "Libint2-export.tgz" "${EXPORT_INSTALL_PREFIX}"
    WORKING_DIRECTORY ${Libint2_BINARY_DIR}
    COMMENT "Building export tarball"
    ALWAYS TRUE
    EXCLUDE_FROM_MAIN FALSE)

ExternalProject_Add(library
    DEPENDS export
    BUILD_ALWAYS 1
    # download generalized to generating source by unpacking source tarball
    URL "${Libint2_BINARY_DIR}/Libint2-export.tgz"
    DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/src/library-prefix/src/library"
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
    TEST_COMMAND ${CMAKE_CTEST_COMMAND})
