include(ExternalProject)
include(GNUInstallDirs)

## correcting for this repo being one off from Psi4 conventions
#math(EXPR MAX_AM_ERI "${MAX_AM_ERI}-1")

ExternalProject_Add(generator
    BUILD_ALWAYS 1
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin/libint
    CMAKE_ARGS -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
               -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DLibint2_BINARY_DIR=${Libint2_BINARY_DIR}
               -DLibint2_ROOT=${Libint2_SOURCE_DIR}
               -DMPFR_ROOT=${MPFR_ROOT}
               -DBOOST_ROOT=${BOOST_ROOT}
    CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS})

ExternalProject_Add(library
    DEPENDS generator
    BUILD_ALWAYS 1
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/libint
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/libint-src
    # download generalized to generating source by running build_libint generator
    DOWNLOAD_COMMAND ${STAGED_INSTALL_PREFIX}/bin/build_libint
    LOG_DOWNLOAD 1
    BINARY_DIR ${CMAKE_BINARY_DIR}/libint-src
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
               -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
               -DLibint2_BINARY_DIR=${Libint2_BINARY_DIR}
               -DLibint2_ROOT=${Libint2_SOURCE_DIR}
               -DLIBINT_MAJOR_SOVERSION=${LIBINT_MAJOR_SOVERSION}
               -DBUILD_SHARED=${BUILD_SHARED}
               -DBUILD_STATIC=${BUILD_STATIC}
               -DMPFR_ROOT=${MPFR_ROOT}
               -DBOOST_ROOT=${BOOST_ROOT}
    CMAKE_CACHE_ARGS -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS})
