cmake_minimum_required(VERSION 3.19)  # string(json
cmake_policy(SET CMP0074 NEW)
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    # Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
    cmake_policy(SET CMP0135 NEW)
endif()
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.25.0")
    # needed by DynamicVersion
    cmake_policy(SET CMP0140 NEW)
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
set(${PROJECT_NAME}_LICENSE "GPL-3.0 for generator; LGPL-3.0 for generated")

# along with project(... VERSION) above scanned by dynamic_version() from `git
#   describe`, these are the authoritative version source (formerly in configure.ac)
set(LIBINT_BUILDID "post999")
set(LIBINT_SOVERSION "2:3:0")
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
   If ERI3 ints are enabled, specifing values for each derivative level as a
   semicolon-separated string also controls the AM of the paired centers." 4)
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
   This option controls only the single fitting center; the paired centers use WITH_MAX_AM." -1)
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

######################## Process & Validate Options ###########################
include(FeatureSummary)
include(int_am)

booleanize01(ERI3_PURE_SH)
booleanize01(ERI2_PURE_SH)
booleanize01(DISABLE_ONEBODY_PROPERTY_DERIVS)
booleanize01(SUPPORT_T1G12)

################################## Main Project #################################

configure_file(include/libint2/config.h.cmake.in include/libint2/config.h @ONLY)

# STRICTLY TEMPORARY FOR DEMONSTRATION PURPOSES
configure_file(src/lib/libint/configuration.cc configuration.cc @ONLY)

message("")
feature_summary(WHAT ENABLED_FEATURES DESCRIPTION "Libint Enabled features:")
feature_summary(WHAT DISABLED_FEATURES DESCRIPTION "Libint Disabled features:")
