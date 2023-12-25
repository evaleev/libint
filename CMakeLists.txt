cmake_minimum_required(VERSION 3.19)  # string(json
cmake_policy(SET CMP0074 NEW)
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    # Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
    cmake_policy(SET CMP0135 NEW)
endif()

############################# Version and Metadata #############################

# can't use PROJECT_SOURCE_DIR etc. before project() call
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
include(DynamicVersion)
dynamic_version(
  PROJECT_PREFIX Libint2Compiler_
  GIT_ARCHIVAL_FILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/.git_archival.txt
  OUTPUT_COMMIT LibintRepository_COMMIT
  OUTPUT_VERSION LibintRepository_VERSION
  OUTPUT_DESCRIBE LibintRepository_DESCRIBE
  OUTPUT_DISTANCE LibintRepository_DISTANCE
  )

project(
  Libint2Compiler
  VERSION ${LibintRepository_VERSION}
  DESCRIPTION
    "A library for the evaluation of molecular integrals of many-body operators over Gaussian functions"
  HOMEPAGE_URL "http://libint.valeyev.net"
  LANGUAGES CXX
  )
  # * http://libint.valeyev.net/ redirects to https://github.com/evaleev/libint

set(${PROJECT_NAME}_AUTHORS "Edward F. Valeev")
set(${PROJECT_NAME}_LICENSE "GPL-3.0-only for generator; LGPL-3.0-only for generated library")

# along with project(... VERSION) above scanned by dynamic_version() from `git
#   describe`, these are the authoritative version source (formerly in configure.ac)
set(LIBINT_BUILDID "post999")
set(LIBINT_SO_VERSION "2:3:0")
set(LIBINT_DOI "10.5281/zenodo.10369117")  # 2.8.0

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
#    >>> cmake --build build --target install

# The Libint build is structured into three parts:
#
# * generator/compiler
#   - (1) build src/bin/libint/ into compiler executable `build_libint`
#   - pretty quick, runs in parallel
#   - consumes all the enable/max/opt/orderings integral options
#   - (2) optionally testable
# * export
#   - (3) run `build_libint` to generate library source (C++) files (that upon
#     compilation can become a Libint2 library) and combine them with other
#     static source files in src/lib/libint/ and general static files (e.g.,
#     include/ and docs/) into an independent tarball ready for distribution
#     (with its own CMake configuration, tests, etc.).
#   - really slow for non-trivial angular momenta; runs in serial
#   - consumes no options
#   - build target `export` to stop after this step and collect source tarball
# * library
#   - can be built as a subproject (FetchContent) or completely insulated (bare
#     ExternalProject; default; -or- a tarball start). For FetchContent, must
#     build libint-library-export target before library build targets appear
#   - (4) unpack the export tarball and build the library and install into \<build\>/library-install-stage/
#   - duration depends on number of integrals requested; runs in parallel
#   - consumes language-interface and the CMAKE_INSTALL_[DATA|INCLUDE|LIB]DIR paths options
#   - the default build target includes this final library build
#   - (5) optionally testable
#   - (6) install into CMAKE_INSTALL_PREFIX
#   - (7) optional Python build alongside library or afterwards. Optional testing requires library install

#################################### Guide #####################################

# See INSTALL.md for elaboration of steps above, options below, & translations from libtool.

################################### Options ####################################
include(options)
include(GNUInstallDirs)
include(CTest)
message(STATUS "Building using CMake ${CMAKE_VERSION} Generator ${CMAKE_GENERATOR}")

#  <<<  General  >>>

### compiler-only

### library-only
option_with_default(LIBINT2_ENABLE_PYTHON
  "Build Python bindings (requires Python and Pybind11 and Eigen3)" OFF)
option_with_default(LIBINT2_PREFIX_PYTHON_INSTALL
  "For LIBINT2_ENABLE_PYTHON=ON, whether to install the Python module in the Linux manner to CMAKE_INSTALL_PREFIX or to not install it. See target libint2-python-wheel for alternate installation in the Python manner to Python_EXECUTABLE's site-packages." OFF)

option_with_print(LIBINT2_ENABLE_MPFR
  "Use GNU MPFR library for high-precision testing (EXPERTS ONLY). Consumed at library build-time." OFF)
## next one defined by `include(CTest)`

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
   These derivatives are disabled by default to save compile time. (enable with OFF)
   Note that the libtool build won't enable this- if forcibly enabled, build_libint balks." ON)
option_with_print(ENABLE_T1G12_SUPPORT
  "Enable Ti,G12 integrals when G12 integrals are enabled. Irrelevant when `ENABLE_G12=OFF`. (disable with OFF)" ON)

#  <<<  Ordering Conventions  >>>

option_with_default(LIBINT2_SHGAUSS_ORDERING
  "Ordering for shells of solid harmonic Gaussians:
    standard -- standard ordering (-l, -l+1 ... l)
    gaussian -- the Gaussian ordering (0, 1, -1, 2, -2, ... l, -l)
   See https://github.com/evaleev/libint/blob/master/INSTALL.md#solid-harmonic-ordering-scope-and-history ." standard)
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

#  <<<  How High Angular Momentum  >>>

#  example for "semicolon-separated string": `-DENABLE_ERI3=2 -DWITH_ERI3_MAX_AM="5;4;3"`

#  special considerations for high-AM library builds:
#  * high MAX_AM generates a large number of source files. If unity builds are disabled, more than
#    ~20k files may require `ulimit -s 65535` for linking the library target on Linux to avert
#    "ld: Argument list too long".
#  * Ninja builds use beyond max threads and can run out of memory, resulting in errorless stops or
#    "CMake Error: Generator: execution of make failed". Throttle it to physical threads with
#    `export CMAKE_BUILD_PARALLEL_LEVEL=N`.

option_with_default(WITH_MAX_AM
  "Support Gaussians of angular momentum up to N.
   Can specify values for each derivative level as a semicolon-separated string.
   If ERI3 ints are enabled, this option also controls the AM of the paired centers." 4)
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
   This option controls only the single fitting center. The paired centers use WITH_MAX_AM." -1)
option_with_default(WITH_ERI3_OPT_AM
  "Optimize 3-center ERIs maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative level as a semicolon-separated string. (default: (max_am/2)+1)" -1)
option_with_print(ERI3_PURE_SH
  "Assume the 'unpaired' center of 3-center ERIs will be transformed to pure solid harmonics" OFF)

option_with_default(WITH_ERI2_MAX_AM
  "Support 2-center ERIs for Gaussians of angular momentum up to N.
    Can specify values for each derivative level as a semicolon-separated string. (default: max_am)" -1)
option_with_default(WITH_ERI2_OPT_AM
  "Optimize 2-center ERIs maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative level as a semicolon-separated string. (default: (max_am/2)+1)" -1)
option_with_print(ERI2_PURE_SH
  "Assume the 2-center ERIs will be transformed to pure solid harmonics" OFF)

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
  "(EXPERT) Compute ERIs using the following strategy. (0 for OS, 1 for HGP, 2 for HL)" 1)
option_with_print(LIBINT_USE_COMPOSITE_EVALUATORS
  "Libint will use composite evaluators (i.e. every evaluator will compute one integral type only)" ON)
option_with_print(LIBINT_SINGLE_EVALTYPE
  "Generate single evaluator type (i.e. all tasks use the same evaluator). OFF is NYI" ON)
option_with_default(LIBINT_ENABLE_UNROLLING
  "Unroll shell sets into integrals (will unroll shell sets larger than N) (0 for never, N for N, 1000000000 for always)" 100)
option_with_default(LIBINT_ALIGN_SIZE
  "(EXPERT) if posix_memalign is available, this will specify alignment of Libint data, in units of
   sizeof(LIBINT2_REALTYPE). Default is to use built-in heuristics: system-determined for vectorization off (default) or veclen * sizeof(LIBINT2_REALTYPE) for vectorization on." 0)
mark_as_advanced(LIBINT_ALIGN_SIZE)
option_with_default(LIBINT2_REALTYPE
  "Specifies the floating-point data type used by the library. Consumed at library build-time." double)
option_with_print(LIBINT_USER_DEFINED_REAL_INCLUDES
  "Additional #includes necessary to use the real type." OFF)
include(int_userreal)
option_with_print(LIBINT_GENERATE_FMA
  "Generate FMA (fused multiply-add) instructions (to benefit must have FMA-capable hardware and compiler)" OFF)
option_with_print(LIBINT_ENABLE_GENERIC_CODE
  "Use manually-written generic code" OFF)
option_with_print(LIBINT_API_PREFIX
  "Prepend this string to every name in the library API (except for the types)." OFF)
option_with_print(LIBINT_VECTOR_LENGTH
  "Compute integrals in vectors of length N." OFF)
option_with_default(LIBINT_VECTOR_METHOD
  "Specifies how to vectorize integrals. Irrelevant when `LIBINT_VECTOR_LENGTH=OFF. Allowed values are 'block' (default), and 'line'." block)
option_with_print(LIBINT_ACCUM_INTS
  "Accumulate integrals to the buffer, rather than copy (OFF for copy, ON for accum)." OFF)
option_with_print(LIBINT_FLOP_COUNT
  "Support (approximate) FLOP counting by the library. (Generated code will require C++11!)" OFF)
option_with_print(LIBINT_PROFILE
  "Turn on profiling instrumentation of the library. (Generated code will require C++11!)" OFF)



######################## Process & Validate Options ###########################
include(CheckFunctionExists)
include(CheckIncludeFileCXX)
include(FeatureSummary)
include(int_orderings)
include(int_am)

check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
if (NOT HAVE_POSIX_MEMALIGN)
    message(FATAL_ERROR "did not find posix_memalign ... this SHOULD NOT happen. Cannot proceed.")
endif()

check_include_file_cxx(stdint.h HAVE_STDINT_H)  # limits.h?

if (cxx_std_11 IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(LIBINT_HAS_CXX11 1)
endif()

booleanize01(ERI3_PURE_SH)
booleanize01(ERI2_PURE_SH)
booleanize01(LIBINT_SINGLE_EVALTYPE)
booleanize01(LIBINT_CONTRACTED_INTS)
booleanize01(DISABLE_ONEBODY_PROPERTY_DERIVS)
booleanize01(LIBINT_GENERATE_FMA)
booleanize01(LIBINT_ENABLE_GENERIC_CODE)
booleanize01(SUPPORT_T1G12)
booleanize01(LIBINT_ACCUM_INTS)
booleanize01(LIBINT_FLOP_COUNT)
booleanize01(LIBINT_PROFILE)


################################## Dependencies #################################

# See notes at https://github.com/evaleev/libint/blob/master/INSTALL.md#prerequisites

if (LIBINT2_ENABLE_MPFR)
    # mpfr detected at appropriate time for library (CMakeLists.txt.export), but prechecking here
    find_package(Multiprecision MODULE REQUIRED COMPONENTS gmpxx mpfr)
    set(LIBINT_HAS_MPFR 1)
else()
    find_package(Multiprecision MODULE REQUIRED COMPONENTS gmpxx)
endif()

get_property(_loc TARGET Multiprecision::gmp PROPERTY LOCATION)
message(STATUS "${Cyan}Found GMP${ColourReset}: ${_loc}")
get_property(_loc TARGET Multiprecision::gmpxx PROPERTY LOCATION)
message(STATUS "${Cyan}Found GMPXX${ColourReset}: ${_loc}")
if (TARGET Multiprecision::mpfr)
    get_property(_loc TARGET Multiprecision::mpfr PROPERTY LOCATION)
    message(STATUS "${Cyan}Found MPFR${ColourReset}: ${_loc} (found version ${MPFR_VERSION})")
endif()

find_package(Boost 1.57 REQUIRED)
if (TARGET Boost::headers)
    set(LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS 1)
endif()

# deferring find_package(Eigen3) to library (CMakeLists.txt.export)


################################## Main Project #################################

set(EXPORT_STAGE_DIR ${PROJECT_BINARY_DIR}/libint-${LIBINT_EXT_VERSION})

configure_file(
  cmake/modules/int_computed.cmake.in
  cmake/modules/int_computed.cmake
  @ONLY)

# CMake data transmitted to C++ via config.h for generator/compiler (_EXPORT_MODE=0).
#   Same info is positioned for the library export, but _EXPORT_MODE=1 turns on
#   inclusion of config2.h to be filled in later by user re-set-able options at
#   library build time.
set(_EXPORT_MODE 0)
configure_file(
  include/libint2/config.h.cmake.in
  include/libint2/config.h
  @ONLY)
set(_EXPORT_MODE 1)
configure_file(
  include/libint2/config.h.cmake.in
  ${EXPORT_STAGE_DIR}/include/libint2/config.h
  @ONLY)
configure_file(
  include/libint2/config2.h.cmake.in
  ${EXPORT_STAGE_DIR}/include/libint2/config2.h.cmake.in
  COPYONLY)

add_subdirectory(src)

message("")
feature_summary(WHAT ENABLED_FEATURES DESCRIPTION "Libint Enabled features:")
feature_summary(WHAT DISABLED_FEATURES DESCRIPTION "Libint Disabled features:")
