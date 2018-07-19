include(ExternalProject)
include(GNUInstallDirs)

# correcting for this repo being one off from Psi4 conventions
math(EXPR MAX_AM_ERI "${MAX_AM_ERI}-1")

ExternalProject_Add(libint_compiler
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin/libint
    CMAKE_ARGS -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
               -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DMAX_AM_ERI=${MAX_AM_ERI}
    CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS})

ExternalProject_Add(libint_library
    DEPENDS libint_compiler
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/libint
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/libint-src
    # download generalized to generating source by running libint_compiler
    DOWNLOAD_COMMAND ${STAGED_INSTALL_PREFIX}/bin/libint_compiler
    LOG_DOWNLOAD 1
    BINARY_DIR ${CMAKE_BINARY_DIR}/libint-src
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
               -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
               -DMAX_AM_ERI=${MAX_AM_ERI}
               -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
               -DBUILD_FPIC=${BUILD_FPIC}
               -DLIBC_INTERJECT=${LIBC_INTERJECT}
    CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                     -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS})

ExternalProject_Add(libderiv_compiler
    DEPENDS libint_library
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin/libderiv
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DMAX_AM_ERI=${MAX_AM_ERI}
               -DLibintint_DIR=${STAGED_INSTALL_PREFIX}/share/cmake/Libint
    CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS})

ExternalProject_Add(libderiv_library
    DEPENDS libderiv_compiler
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/libderiv
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/libderiv-src
    # download generalized to generating source by running libderiv_compiler
    DOWNLOAD_COMMAND ${STAGED_INSTALL_PREFIX}/bin/libderiv_compiler
    LOG_DOWNLOAD 1
    BINARY_DIR ${CMAKE_BINARY_DIR}/libderiv-src
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
               -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
               -DMAX_AM_ERI=${MAX_AM_ERI}
               -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
               -DBUILD_FPIC=${BUILD_FPIC}
               -DLIBC_INTERJECT=${LIBC_INTERJECT}
               -DLibintint_DIR=${STAGED_INSTALL_PREFIX}/share/cmake/Libint
               -DMERGE_LIBDERIV_INCLUDEDIR=${MERGE_LIBDERIV_INCLUDEDIR}
    CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                     -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS})

# Note that libint_library and libderiv_library could build in parallel
#   if (1) the code generation step of libint_library was moved to the end
#   of libint_compiler and (2) libderiv_compiler was given a way to detect
#   the generated files. That code generation step produces libint.h and
#   that's what libderiv_compiler is waiting for.
