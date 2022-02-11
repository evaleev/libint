cmake_minimum_required(VERSION 3.18)  # LIST SORT NATURAL
cmake_policy(SET CMP0074 NEW)

############################# Version and Metadata #############################

project(
  Libint2Compiler
  VERSION 2.7.1
  DESCRIPTION "High-performance library for computing Gaussian integrals in quantum mechanics"
  HOMEPAGE_URL "https://github.com/evaleev/libint"
  LANGUAGES CXX
  )

set(${PROJECT_NAME}_AUTHORS "Edward F. Valeev")
set(${PROJECT_NAME}_LICENSE "GPL-3.0 for generator; LGPL-3.0 for generated")

# along with project(... VERSION) above, these are the authoritative version source (formerly in configure.ac)
set(LIBINT_BUILDID "")
set(LIBINT_SOVERSION "2:3:0")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
include(int_versions)

set(L2 Libint2)  # Namespace
set(pnv libint2) # projectnameversion

################################### Overview ###################################

# CMake build overview:
#
#    >>> ls
#    cmake/  COPYING  src/  tests/  ...
#    >>> cmake -S. -Bbuild -GNinja -DCMAKE_INSTALL_PREFIX=/path/to/install-libint ...
#    ...
#    -- Generating done
#    -- Build files have been written to: /current/dir/build
#    >>> cmake --build build --target install -j`getconf _NPROCESSORS_ONLN`

# The build is structured into three parts:
#
# * generator
#   - build src/bin/libint/ into executable `build_libint`
#   - pretty quick, runs in parallel
#   - consumes all the enable/max/opt integral options and all orderings options except solid harmonic
# * export
#   - run `build_libint` to generate library source (C++) files (that upon
#     compilation can become a Libint2 library) and combine them with other
#     static source files in src/lib/libint/ and general static files (e.g.,
#     include/ and docs/) into an independent tarball ready for distribution.
#   - really slow for non-trivial angular momenta, runs in serial
#   - consumes no options
#   - build target `export` to stop after this step and collect source tarball
# * library
#   - can be build as a subproject (FetchContent) or completely insulated (bare ExternalProject; default)
#   - if building via bare ExternalProject:
#     - unpack the export tarball and build the library and install into stage/
#     - duration depends on number of integrals requested, runs in parallel
#     - consumes solid harmonic ordering and the CMAKE_INSTALL_[DATA|INCLUDE|LIB]DIR
#     - the default target build includes this final library build
#   - if building via FetchContent:
#     - must build libint-library-export target before library *build* targets appear

################################### Options ####################################
include(options)
include(GNUInstallDirs)
include(CTest)
message(STATUS "Building using CMake ${CMAKE_VERSION}")
message(STATUS "Generating using ${CMAKE_GENERATOR}: ${CMAKE_GENERATOR_INSTANCE}")

#  <<<  General  >>>

option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)

### compiler-only
option_with_default(LIBINT_BUILD_LIBRARY_AS_SUBPROJECT "[EXPERT] Build generated library as a subproject; if FALSE will configure and build separately" OFF)

### library-only
option_with_default(REQUIRE_CXX_API "C++11 Libint API: define library targets + test (requires Eigen3; Boost is optional but strongly recommended)" ON)
option_with_default(ENABLE_FORTRAN "Build Fortran03+ Libint interface (requires Fortran)" OFF)

option_with_print(BUILD_SHARED_LIBS "Build Libint library as shared, not static" OFF)
option_with_print(LIBINT2_BUILD_SHARED_AND_STATIC_LIBS "Build both shared and static Libint libraries in one shot. Uses -fPIC." OFF)

#  <<<  Which Integrals Classes, Which Derivative Levels  >>>

option_with_default(ENABLE_ONEBODY
  "Compile with support for up to N-th derivatives of 1-body integrals (-1 for OFF)" 0)
option_with_default(ENABLE_ERI
  "Compile with support for up to N-th derivatives of 4-center electron repulsion integrals (-1 for OFF)" 0)
option_with_default(ENABLE_ERI3
  "Compile with support for up to N-th derivatives of 3-center electron repulsion integrals (-1 for OFF)" -1)
option_with_default(ENABLE_ERI2
  "Compile with support for up to N-th derivatives of 2-center electron repulsion integrals (-1 for OFF)" -1)
option_with_default(ENABLE_G12
  "Compile with support for N-th derivatives of MP2-F12 energies with Gaussian factors (-1 for OFF)" -1)
option_with_default(ENABLE_G12DKH
  "Compile with support for N-th derivatives of DKH-MP2-F12 energies with Gaussian factors (-1 for OFF)" -1)

option_with_print(DISABLE_ONEBODY_PROPERTY_DERIVS
  "Disable geometric derivatives of 1-body property integrals (all but overlap, kinetic, elecpot).
   These derivatives are disabled by default to save compile time. (enable with OFF)" ON)
option_with_print(ENABLE_T1G12_SUPPORT
  "Enable Ti,G12 integrals when G12 integrals are enabled. Irrelevant when `ENABLE_G12=OFF`. (disable with OFF)" ON)

#  <<<  Ordering Conventions  >>>

option_with_default(LIBINT2_SHGAUSS_ORDERING
  "Ordering for shells of solid harmonic Gaussians:
    standard -- standard ordering (-l, -l+1 ... l)
    gaussian -- the Gaussian ordering (0, 1, -1, 2, -2, ... l, -l)" standard)
option_with_default(LIBINT2_CARTGAUSS_ORDERING
  "Orderings for shells of cartesian Gaussians:
    standard -- standard ordering (xxx, xxy, xxz, xyy, xyz, xzz, yyy, ...)
    intv3  -- intv3 ordering (yyy, yyz, yzz, zzz, xyy, xyz, xzz, xxy, xxz, xxx)
    gamess -- GAMESS ordering (xxx, yyy, zzz, xxy, xxz, yyx, yyz, zzx, zzy, xyz)
    orca -- ORCA ordering (hydrid between GAMESS and standard)
    bagel -- axis-permuted version of intv3 (xxx, xxy, xyy, yyy, xxz, xyz, yyz, xzz, yzz, zzz)" standard)
option_with_default(LIBINT2_SHELL_SET
  "Support computation of shell sets sets subject to these restrictions:
    standard -- standard ordering:
      for (ab|cd):
        l(a) >= l(b),
        l(c) >= l(d),
        l(a)+l(b) <= l(c)+l(d)
      for (b|cd):
        l(c) >= l(d)
    orca -- ORCA ordering:
      for (ab|cd):
        l(a) <= l(b),
        l(c) <= l(d),
        l(a) < l(c) || (l(a) == l(c) && l(b) < l(d))
      for (b|cd):
        l(c) <= l(d)" standard)
option_with_print(ERI3_PURE_SH
  "Assume the 'unpaired' center of 3-center ERIs will be transformed to pure solid harmonics" OFF)
option_with_print(ERI2_PURE_SH
  "Assume the 2-center ERIs will be transformed to pure solid harmonics" OFF)

#  <<<  How High Angular Momentum  >>>

#  example for "semicolon-separated string": `-DENABLE_ERI3=2 -DWITH_ERI3_MAX_AM="5;4;3"`

#  high MAX_AM generating >20k files may require `ulimit -s 65535` for linking library target on Linux to avert "ld: Argument list too long"

option_with_default(WITH_MAX_AM
  "Support Gaussians of angular momentum up to N.
   Can specify values for each derivative level as a semicolon-separated string" 4)
option_with_default(WITH_OPT_AM
  "Optimize maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative level as a semicolon-separated string. (default: (libint_max_am/2)+1)" -1)

option_with_default(MULTIPOLE_MAX_ORDER
  "Maximum order of spherical multipole integrals. There is no maximum" 4)
option_with_default(WITH_ONEBODY_MAX_AM
  "Support 1-body ints for Gaussians of angular momentum up to N.
   Can specify values for each derivative level as a semicolon-separated string. (default: max_am)" -1)
option_with_default(WITH_ONEBODY_OPT_AM
  "Optimize 1-body ints maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative level as a semicolon-separated string (default: (max_am/2)+1)" -1)

option_with_default(WITH_ERI_MAX_AM
  "Support 4-center ERIs for Gaussians of angular momentum up to N.
   Can specify values for each derivative level as a semicolon-separated string. (default: max_am)" -1)
option_with_default(WITH_ERI_OPT_AM
  "Optimize 4-center ERIs maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative level as a semicolon-separated string (default: (max_am/2)+1)" -1)

option_with_default(WITH_ERI3_MAX_AM
  "Support 3-center ERIs for Gaussians of angular momentum up to N.
   Can specify values for each derivative level as a semicolon-separated string. (default: max_am)
   Set WITH_MAX_AM to mirror this option." -1)
option_with_default(WITH_ERI3_OPT_AM
  "Optimize 3-center ERIs maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative level as a semicolon-separated string. (default: (max_am/2)+1)" -1)

option_with_default(WITH_ERI2_MAX_AM
  "Support 2-center ERIs for Gaussians of angular momentum up to N.
    Can specify values for each derivative level as a semicolon-separated string. (default: max_am)" -1)
option_with_default(WITH_ERI2_OPT_AM
  "Optimize 2-center ERIs maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative level as a semicolon-separated string. (default: (max_am/2)+1)" -1)

option_with_default(WITH_G12_MAX_AM
  "Support integrals for G12 methods of angular momentum up to N. (default: max_am)" -1)
option_with_default(WITH_G12_OPT_AM
  "Optimize G12 integrals for up to angular momentum N (N <= max-am). (default: (max_am/2)+1)" -1)

option_with_default(WITH_G12DKH_MAX_AM
  "Support integrals for relativistic G12 methods of angular momentum up to N. (default: max_am)" -1)
option_with_default(WITH_G12DKH_OPT_AM
  "Optimize G12DKH integrals for up to angular momentum N (N <= max-am). (default: (max_am/2)+1)" -1)

#  <<<  Miscellaneous  >>>

option_with_print(LIBINT_CONTRACTED_INTS
  "Turn on support for contracted integrals." ON)
option_with_default(LIBINT_ERI_STRATEGY
  "Compute ERIs using the following strategy. This option is for experts ONLY." 1)
option_with_print(LIBINT_USE_COMPOSITE_EVALUATORS
  "Libint will use composite evaluators (i.e. every evaluator will compute one integral type only)" ON)
option_with_print(LIBINT_SINGLE_EVALTYPE
  "Generate single evaluator type (i.e. all tasks use the same evaluator)" ON)
option_with_default(LIBINT_ENABLE_UNROLLING
  "Unroll shell sets into integrals (will unroll shell sets larger than N) (no->0, yes->1000000000)" 100)
option_with_default(LIBINT_ALIGN_SIZE
  "(EXPERT) if posix_memalign is available, this will specify alignment of Libint data, in units of
   sizeof(LIBINT2_REALTYPE). Default is to use built-in heuristics. LAB: LIBINT for bin, LIBINT2 for lib?" 0)
option_with_print(LIBINT_GENERATE_FMA
  "Generate FMA (fused multiply-add) instructions (to benefit must have FMA-capable hardware and compiler)" OFF)
option_with_print(LIBINT_ENABLE_GENERIC_CODE
  "Use manually-written generic code" OFF)
#option_with_flags(ENABLE_XHOST "Enables processor-specific optimization (with MSVC, it enables AVX2 instructions)" ON
#                  "-xHost" "-march=native" "/arch:AVX2")
option_with_default(BUILD_FPIC "Compile lone static library with position independent code (all libraries are fpic when BUILD_SHARED=ON)" OFF)
option_with_print(ENABLE_MPFR "Use GNU MPFR library for high-precision testing (EXPERTS ONLY)" OFF)
# next one defined by `include(CTest)`
message(STATUS "Showing option BUILD_TESTING: ${BUILD_TESTING}")

#  <<<  Path  >>>

# * only relevant to library
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

#CONFIGURE    --enable-accum-ints     Accumulate integrals to the buffer, rather than
#CONFIGURE                            copy.
#CONFIGURE    --enable-flop-counter   Support (approximate) FLOP counting by the library.
#CONFIGURE                            CXXGEN must support C++11!
#CONFIGURE    --enable-profile        Turn on profiling instrumentation of the library.
#CONFIGURE                            CXXGEN must support C++11!
#CONFIGURE  
#CONFIGURE  Optional Packages:
#CONFIGURE    --with-api-prefix       Prepend this string to every name in the library API
#CONFIGURE                            (except for the types).
#CONFIGURE    --with-vector-length=N  Compute integrals in vectors of length N.
#CONFIGURE    --with-vector-method    Specifies how to vectorize integrals. Allowed values
#CONFIGURE                            are 'block' (default), and 'line'.
#CONFIGURE    --with-build-id         Gives an identifier for the build.
#CONFIGURE    --with-libint-exportdir Specifies export directory name. Default is
#CONFIGURE                            libint-<version>.
#CONFIGURE    --with-pkgconfigdir     Specifies the location to put pkg-config's data
#CONFIGURE                            file. Default is $(libdir)/pkgconfig.
#CONFIGURE    --with-real-type        Compute all integrals using the given real type           LIBINT_USER_DEFINED_REAL
#CONFIGURE                            (default: double).
#CONFIGURE    --with-real-type-includes                                                         LIBINT_USER_DEFINED_REAL_INCLUDES
#CONFIGURE                            Additional includes necessary to use the real type
#CONFIGURE                            (default: none).

######################## Process & Validate Options ###########################
include(autocmake_safeguards)
include(CheckFunctionExists)
include(int_orderings)
include(int_am)

check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
if (HAVE_POSIX_MEMALIGN)
    set(LIBINT2_ALIGN_SIZE "0" CACHE STRING "(EXPERT) if posix_memalign is available, this will specify alignment of Libint data, in units of sizeof(LIBINT2_REALTYPE). Default is to use built-in heuristics")
endif()

# C++11 required for now. enforced in targets
set(LIBINT_HAS_CXX11 1)
set(LIBINT_USE_TYPE_TRAITS 1)
set(LIBINT_USE_MEMORY 1)
set(LIBINT_USE_ARRAY 1)
set(LIBINT_HAS_STD_TYPE_TRAITS 1)
set(LIBINT_HAS_STD_SHARED_PTR 1)
set(LIBINT_HAS_STD_MAKE_SHARED 1)
set(LIBINT_HAS_STD_ARRAY 1)
set(LIBINT_ARRAY_HAS_FILL 1)

find_package(Boost 1.57 REQUIRED)
if (TARGET Boost::headers)
    set(LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS 1)
endif()

unset(_temp_h CACHE)
find_path(_temp_h NAMES stdint.h)
if(_temp_h)
    set(HAVE_STDINT_H 1)
endif()

booleanize01(ERI3_PURE_SH)
booleanize01(ERI2_PURE_SH)
booleanize01(LIBINT_SINGLE_EVALTYPE)
booleanize01(LIBINT_USE_COMPOSITE_EVALUATORS)
booleanize01(LIBINT_CONTRACTED_INTS)
booleanize01(DISABLE_ONEBODY_PROPERTY_DERIVS)
booleanize01(LIBINT_GENERATE_FMA)
booleanize01(LIBINT_ENABLE_GENERIC_CODE)
booleanize01(SUPPORT_T1G12)

################################# Main Project #################################
set(EXPORT_STAGE_DIR ${PROJECT_BINARY_DIR}/libint-${LIBINT_EXT_VERSION})
set(STAGED_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/library-install-stage")

configure_file(cmake/modules/int_computed.cmake.in cmake/modules/int_computed.cmake @ONLY)
set(_EXPORT_MODE 0)
configure_file(include/libint2/config.h.cmake.in include/libint2/config.h @ONLY)
set(_EXPORT_MODE 1)
configure_file(include/libint2/config.h.cmake.in ${EXPORT_STAGE_DIR}/include/libint2/config.h @ONLY)
configure_file(include/libint2/config2.h.cmake.in ${EXPORT_STAGE_DIR}/include/libint2/config2.h.cmake.in COPYONLY)

add_subdirectory(src)

add_subdirectory(doc)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
