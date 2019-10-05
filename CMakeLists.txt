cmake_minimum_required(VERSION 3.8)  #TODO

project(Libint2
        LANGUAGES CXX)
set(Libint2_AUTHORS      "Edward F. Valeev")
set(Libint2_DESCRIPTION  "High-performance library for computing Gaussian integrals in quantum mechanics")
set(Libint2_URL          "https://github.com/evaleev/libint")
set(Libint2_LICENSE      "GPL-3.0 for generator; LGPL-3.0 for generated")

cmake_policy(SET CMP0074 NEW)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake
                              ${PROJECT_SOURCE_DIR}/export/cmake/modules)

################################### Options ####################################
include(int_versions)
include(psi4OptionsTools)

#  <<< General  >>>

option_with_default(CMAKE_BUILD_TYPE "Build type" Release)
option_with_print(BUILD_SHARED "Build at least shared library" OFF)
option_with_print(BUILD_STATIC "Build at least static library" ON)

#  <<<  Which Integrals Classes, Which Derivative Levels  >>>

option_with_default(ENABLE_ONEBODY
  "Compile with support for up to N-th derivatives of 1-body integrals (-1 for OFF)" 0)
option_with_default(ENABLE_ERI
  "Compile with support for up to N-th derivatives of electron repulsion integrals (-1 for OFF)" 0)
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

# example for "semicolon-separated string": `-DENABLE_ERI3=2 -DWITH_ERI3_MAX_AM="5;4;3"`

option_with_default(WITH_MAX_AM
  "Support Gaussians of angular momentum up to N.
   Can specify values for each derivative levels as a semicolon-separated string" 4)
option_with_default(WITH_OPT_AM
  "Optimize maximally for up to angular momentum N (N <= max-am). (default: (libint_max_am/2)+1)" -1)

option_with_default(MULTIPOLE_MAX_ORDER
  "Maximum order of spherical multipole integrals. There is no maximum" 4)
option_with_default(WITH_ONEBODY_MAX_AM
  "Support 1-e ints for Gaussians of angular momentum up to N.
   Can specify values for each derivative as a semicolon-separated string. (default: (libint_max_am/2)+1)" -1)
option_with_default(WITH_ONEBODY_OPT_AM
  "Optimize 1-e ints maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative levels as a semicolon-separated string" -1)

option_with_default(WITH_ERI_MAX_AM
  "Support ERIs for Gaussians of angular momentum up to N.
   Can specify values for each derivative levels as a semicolon-separated string" -1)
option_with_default(WITH_ERI_OPT_AM
  "Optimize ERIs maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative levels as a semicolon-separated string (default: (max_am/2)+1)" -1)

option_with_default(WITH_ERI3_MAX_AM
  "Support 3-center ERIs for Gaussians of angular momentum up to N.
   Can specify values for each derivative levels as a semicolon-separated string" -1)
option_with_default(WITH_ERI3_OPT_AM
  "Optimize 3-center ERIs maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative levels as a semicolon-separated string (default: (max_am/2)+1)" -1)

option_with_default(WITH_ERI2_MAX_AM
  "Support 2-center ERIs for Gaussians of angular momentum up to N.
    Can specify values for each derivative levels as a semicolon-separated string" -1)
option_with_default(WITH_ERI2_OPT_AM
  "Optimize 2-center ERIs maximally for up to angular momentum N (N <= max-am).
   Can specify values for each derivative levels as a semicolon-separated string" -1)

option_with_default(WITH_G12_MAX_AM
  "Support integrals for G12 methods of angular momentum up to N." -1)
option_with_default(WITH_G12_OPT_AM
  "Optimize G12 integrals for up to angular momentum N (N <= max-am)." -1)

option_with_default(WITH_G12DKH_MAX_AM
  "Support integrals for relativistic G12 methods of angular momentum up to N." -1)
option_with_default(WITH_G12DKH_OPT_AM
  "Optimize G12DKH integrals for up to angular momentum N (N <= max-am)." -1)

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


#option_with_print(BUILD_SHARED_LIBS "Build final library as shared, not static" OFF)
#option_with_default(BUILD_FPIC "Libraries will be compiled with position independent code" ON)
#if(${BUILD_SHARED_LIBS} AND NOT ${BUILD_FPIC})
#    message(FATAL_ERROR "BUILD_SHARED_LIBS ON and BUILD_FPIC OFF are incompatible, as shared library requires position independent code")
#endif()
#option_with_print(ENABLE_GENERIC "Enable mostly static linking in shared library" OFF)
#option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" ON
#                    "-xHost" "-march=native")
#option_with_print(MERGE_LIBDERIV_INCLUDEDIR "Install libderiv headers to libint namespace. Psi4 1.1-1.3 wants ON" OFF)
option_with_print(ENABLE_FORTRAN "Build Fortran 90 interface" OFF)

#CONFIGURE  System types:
#CONFIGURE    --build=BUILD     configure for building on BUILD [guessed]
#CONFIGURE    --host=HOST       cross-compile to build programs to run on HOST [BUILD]
#CONFIGURE    --target=TARGET   configure for building compilers for TARGET [HOST]
#CONFIGURE  
#CONFIGURE    --disable-t1g12-support Libint will not support [Ti,G12] integrals
#CONFIGURE    --enable-generic-code   Use manually-written generic code.
#CONFIGURE    --enable-fma            Generate FMA (fused multiply-add) instructions
#CONFIGURE                            (CXXGEN must support C++11; to benefit must have
#CONFIGURE                            FMA-capable hardware and compiler)
#CONFIGURE    --enable-accum-ints     Accumulate integrals to the buffer, rather than
#CONFIGURE                            copy.
#CONFIGURE    --enable-flop-counter   Support (approximate) FLOP counting by the library.
#CONFIGURE                            CXXGEN must support C++11!
#CONFIGURE    --enable-profile        Turn on profiling instrumentation of the library.
#CONFIGURE                            CXXGEN must support C++11!
#CONFIGURE    --enable-mpfr           Use MPFR library to test Libint integrals in high
#CONFIGURE                            precision.
#CONFIGURE    --enable-shared[=PKGS]  build shared libraries [default=no]
#CONFIGURE    --enable-static[=PKGS]  build static libraries [default=yes]
#CONFIGURE    --enable-fast-install[=PKGS]
#CONFIGURE                            optimize for fast installation [default=yes]
#CONFIGURE    --disable-libtool-lock  avoid locking (might break parallel builds)
#CONFIGURE  
#CONFIGURE  Optional Packages:
#CONFIGURE    --with-PACKAGE[=ARG]    use PACKAGE [ARG=yes]
#CONFIGURE    --without-PACKAGE       do not use PACKAGE (same as --with-PACKAGE=no)
#CONFIGURE    --with-api-prefix       Prepend this string to every name in the library API
#CONFIGURE                            (except for the types).
#CONFIGURE    --with-vector-length=N  Compute integrals in vectors of length N.
#CONFIGURE    --with-vector-method    Specifies how to vectorize integrals. Allowed values
#CONFIGURE                            are 'block' (default), and 'line'.
#CONFIGURE    --with-build-id         Gives an identifier for the build.
#CONFIGURE    --with-libint-includedir
#CONFIGURE                            Specifies include file install subdir.
#CONFIGURE    --with-libint-exportdir Specifies export directory name. Default is
#CONFIGURE                            libint-<version>.
#CONFIGURE    --with-incdirs          Specifies include directories (-Idir1 -Idir2).
#CONFIGURE    --with-libs             Specifies libraries (-llib1 -llib2).
#CONFIGURE    --with-libdirs          Specifies library directories (-Ldir1 -Ldir2).
#CONFIGURE    --with-scratchdir       Specifies the location for the machine-generated
#CONFIGURE                            library source and object files.
#CONFIGURE    --with-pkgconfigdir     Specifies the location to put pkg-config's data
#CONFIGURE                            file. Default is $(libdir)/pkgconfig.
#CONFIGURE    --with-cmakedir         Specifies the location to put FindLibint2.cmake
#CONFIGURE                            file. Default is $(libdir)/cmake/libint2.
#CONFIGURE    --with-boost[=ARG]      use Boost library from a standard location
#CONFIGURE                            (ARG=yes), from the specified location (ARG=<path>),
#CONFIGURE                            or disable it (ARG=no) [ARG=yes]
#CONFIGURE    --with-boost-libdir=LIB_DIR
#CONFIGURE                            Force given directory for boost libraries. Note that
#CONFIGURE                            this will override library path detection, so use
#CONFIGURE                            this parameter only if default library detection
#CONFIGURE                            fails and you know exactly where your boost
#CONFIGURE                            libraries are located.
#CONFIGURE    --with-real-type        Compute all integrals using the given real type
#CONFIGURE                            (default: double).
#CONFIGURE    --with-real-type-includes
#CONFIGURE                            Additional includes necessary to use the real type
#CONFIGURE                            (default: none).
#CONFIGURE    --with-pic[=PKGS]       try to use only PIC/non-PIC objects [default=use
#CONFIGURE                            both]
#CONFIGURE    --with-aix-soname=aix|svr4|both
#CONFIGURE                            shared library versioning (aka "SONAME") variant to
#CONFIGURE                            provide on AIX, [default=aix].
#CONFIGURE    --with-gnu-ld           assume the C compiler uses GNU ld [default=no]
#CONFIGURE    --with-sysroot[=DIR]    Search for dependent libraries within DIR (or the
#CONFIGURE                            compiler's sysroot if not specified).
#CONFIGURE  
#CONFIGURE  Some influential environment variables:
#CONFIGURE    CXX         C++ compiler command
#CONFIGURE    CXXFLAGS    C++ compiler flags
#CONFIGURE    LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
#CONFIGURE                nonstandard directory <lib dir>
#CONFIGURE    LIBS        libraries to pass to the linker, e.g. -l<library>
#CONFIGURE    CPPFLAGS    (Objective) C/C++ preprocessor flags, e.g. -I<include dir> if
#CONFIGURE                you have headers in a nonstandard directory <include dir>
#CONFIGURE    CXXCPP      C++ preprocessor
#CONFIGURE    CC          C compiler command
#CONFIGURE    CFLAGS      C compiler flags
#CONFIGURE    LT_SYS_LIBRARY_PATH
#CONFIGURE                User-defined run-time library search path.

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

# req'd for now (no FindGMP.cmake, though one at Elemental project)
find_package(MPFR REQUIRED)
if (TARGET MPFR::MPFRXX)
    set(LIBINT_HAS_MPFR 1)
endif()

find_package(Boost 1.57 REQUIRED)
if (TARGET Boost::boost)  # the header-only target
    set(LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS 1)
endif()

find_package(Eigen3 REQUIRED)
if (TARGET Eigen3::Eigen)
    set(LIBINT_HAS_EIGEN 1)
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

################################# Main Project #################################
set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)

add_subdirectory(src)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(PN ${PROJECT_NAME})

if (ENABLE_FORTRAN)
    enable_language(Fortran)
endif()

configure_file(include/libint2/config.h.cmake.in include/libint2/config.h @ONLY)
set(DATADIR_ABSOLUTE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION})
configure_file(include/libint2/basis.h.in include/libint2/basis.h @ONLY)
##configure_file(cmake/libint.pc.in libint.pc @ONLY)

# <<<  Install  >>>

install(DIRECTORY ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
                  #${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
                  ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
                  ${STAGED_INSTALL_PREFIX}/share
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS)

# <<<  Export Config  >>>

# GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
message("share ${CMAKECONFIG_INSTALL_DIR}")
configure_package_config_file(cmake/${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
                                 VERSION ${${PN}_VERSION}
                                 COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

#install(FILES ${PROJECT_BINARY_DIR}/libint.pc
#        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig/)
