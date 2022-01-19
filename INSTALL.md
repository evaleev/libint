# Synopsis

- configure: `cmake -S /path/to/compiler/source/tree -B build [-Dvar1=value1] [-Dvar2=value2] ...` where
  the optional CMake cache variables `var1`, `var2` are described below. Replace `build` with desired path to
  the build directory, if desired; it will be created, if needed.
- build: `cmake --build build`
- test (optional): `cmake --build build --target check`
- install the library: `cmake --build build --target install`

# Build instructions

For (more) complete instructions please refer to https://github.com/evaleev/libint/wiki

# Appendix


####################

target                      includes steps (see below)
------                      --------------------------
build_libint                1    -    -    -    -    -
check-libint2compiler       1    2    -    -    -    -
export                      1    -    3    -    -    -
library (default)           1    -    3    4    -    -
check                       1    2    3    4    5    -
install                     1    -    3    4    -    6


The build is structured into three parts:

* generator
  - (1) build src/bin/libint/ into generator/compiler executable `build_libint`
  - pretty quick, runs in parallel
  - consumes all the enable/max/opt integral options and all orderings options except solid harmonic
  - (2) optionally testable
* export
  - (3) run `build_libint` to generate library source (C++) files (that upon
    compilation can become a Libint2 library) and combine them with other
    static source files in src/lib/libint/ and general static files (e.g.,
    include/ and docs/) into an independent tarball ready for distribution.
  - really slow for non-trivial angular momenta, runs in serial
  - consumes no options
  - build target `export` to stop after this step and collect source tarball
* library
  - can be build as a subproject (FetchContent) or completely insulated (bare ExternalProject; default)
  - if building via bare ExternalProject:
    - (4) unpack the export tarball and build the library and install into <build>/library-install-stage/
    - duration depends on number of integrals requested, runs in parallel
    - consumes solid harmonic ordering and the CMAKE_INSTALL_[DATA|INCLUDE|LIB]DIR
    - the default build target includes this final library build
    - (5) optionally testable
    - (6) install into CMAKE_INSTALL_PREFIX
  - if building via FetchContent:
    - must build libint-library-export target before library *build* targets appear

-----------------------------------------------------------------------------
## configuring libint compiler

These are the most useful configure options:



###  Which Integrals Classes, Which Derivative Levels

* `ENABLE_ONEBODY` — Compile with support for up to N-th derivatives of 1-body integrals. Use -1 for OFF. [Default=0]
* `ENABLE_ERI` — Compile with support for up to N-th derivatives of 4-center electron repulsion integrals. Use -1 for OFF. [Default=0]
* `ENABLE_ERI3` — Compile with support for up to N-th derivatives of 3-center electron repulsion integrals. Use -1 for OFF. [Default=-1]
* `ENABLE_ERI2` — Compile with support for up to N-th derivatives of 2-center electron repulsion integrals. Use -1 for OFF. [Default=-1]
* `ENABLE_G12` — Compile with support for N-th derivatives of MP2-F12 energies with Gaussian factors. Use -1 for OFF. [Default=-1]
* `ENABLE_G12DKH` — Compile with support for N-th derivatives of DKH-MP2-F12 energies with Gaussian factors. Use -1 for OFF. [Default=-1]

* `DISABLE_ONEBODY_PROPERTY_DERIVS` — Disable geometric derivatives of 1-body property integrals (all but overlap, kinetic, elecpot).
   These derivatives are disabled by default to save compile time. Use OFF to enable. [Default=ON]

###  Ordering Conventions

* `LIBINT2_SHGAUSS_ORDERING` — Ordering for shells of solid harmonic Gaussians. Consumed at library build-time. [Default=standard]
  * `standard` — standard ordering (-l, -l+1 ... l)
  * `gaussian` — the Gaussian ordering (0, 1, -1, 2, -2, ... l, -l)
* `LIBINT2_CARTGAUSS_ORDERING` — Orderings for shells of cartesian Gaussians. Consumed at generator build-time. [Default=standard]
  * `standard` — standard ordering (xxx, xxy, xxz, xyy, xyz, xzz, yyy, ...)
  * `intv3` — intv3 ordering (yyy, yyz, yzz, zzz, xyy, xyz, xzz, xxy, xxz, xxx)
  * `gamess` — GAMESS ordering (xxx, yyy, zzz, xxy, xxz, yyx, yyz, zzx, zzy, xyz)
  * `orca` — ORCA ordering (hydrid between GAMESS and standard)
  * `bagel` — axis-permuted version of intv3 (xxx, xxy, xyy, yyy, xxz, xyz, yyz, xzz, yzz, zzz)
* `LIBINT2_SHELL_SET` — Support computation of shell sets sets subject to these restrictions. Consumed at generator build-time. [Default=standard]
  * `standard` — standard ordering:
      for (ab|cd):
        l(a) >= l(b),
        l(c) >= l(d),
        l(a)+l(b) <= l(c)+l(d)
      for (b|cd):
        l(c) >= l(d)
  * `orca` — ORCA ordering:
      for (ab|cd):
        l(a) <= l(b),
        l(c) <= l(d),
        l(a) < l(c) || (l(a) == l(c) && l(b) < l(d))
      for (b|cd):
        l(c) <= l(d)
* `ERI3_PURE_SH` — Assume the 'unpaired' center of 3-center ERIs will be transformed to pure solid harmonics. [Default=OFF]
* `ERI2_PURE_SH` — Assume the 2-center ERIs will be transformed to pure solid harmonics. [Default=OFF]

###  How High Angular Momentum

* Notes
  * example for "semicolon-separated string": `-DENABLE_ERI3=2 -DWITH_ERI3_MAX_AM="5;4;3"`. cmake configuration prints:

    ```
    -- Setting option ENABLE_ERI3: 2
    -- Setting option WITH_ERI3_MAX_AM: 5;4;3
    ```

  * high MAX_AM generating >20k files may require `ulimit -s 65535` for linking library target on Linux to avert "ld: Argument list too long"

* `WITH_MAX_AM` — Support Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=4]
* `WITH_OPT_AM` — Optimize maximally for up to angular momentum N (N <= WITH_MAX_AM). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `(WITH_MAX_AM/2)+1`]

* `WITH_ONEBODY_MAX_AM` — Support 1-body ints for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ONEBODY_OPT_AM` — Optimize 1-body ints maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI_MAX_AM` — Support 4-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI_OPT_AM` — Optimize 4-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI3_MAX_AM` — Support 3-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI3_OPT_AM` — Optimize 3-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI2_MAX_AM` — Support 2-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI2_OPT_AM` — Optimize 2-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]



### Autotools Update Guide

* Notes
  * When three option names present, they are for libtool+cmake --> cmake+cmake (c. #148 for Psi4 c. 2020-2021) --> cmake+cmake

* `--enable-1body=N` --> `-D ENABLE_ONEBODY`
* `--enable-eri=N` --> `-D ENABLE_ERI=N`
* `--disable-eri` --> `-D ENABLE_ERI=-1`
* `--enable-eri3=N` --> `-D ENABLE_ERI3=N`
* `--enable-eri2=N` --> `-D ENABLE_ERI2=N`

* `--with-shgauss-ordering=label` --> `-D LIBINT2_SHGAUSS_ORDERING=label`
* `--with-cartgauss-ordering=label` --> `-D LIBINT2_CARTGAUSS_ORDERING=label`
* `--with-shell-set=label` --> `-D LIBINT2_SHELL_SET=label`
* `--enable-eri3-pure-sh` --> `-D ERI3_PURE_SH=ON`
* `--enable-eri2-pure-sh` --> `-D ERI2_PURE_SH=ON`

* `--with-max-am=N` --> `-D WITH_MAX_AM=N`
* `--with-max-am=N0,N1,N2` --> `-D WITH_MAX_AM="N0;N1;N2"` (notice semicolons and quotes. This is standard CMake list syntax)
* `--with-opt-am=N` --> `-D WITH_OPT_AM=N`
* `--with-opt-am=N0,N1,N2` --> `-D WITH_OPT_AM="N0;N1;N2"`

* `--with-1body-max-am=N` --> `-D WITH_ONEBODY_MAX_AM=N`
* `--with-1body-max-am=N0,N1,N2` --> `-D WITH_ONEBODY_MAX_AM="N0;N1;N2"`
* `--with-1body-opt-am=N` --> `-D WITH_ONEBODY_OPT_AM=N`
* `--with-1body-opt-am=N0,N1,N2` --> `-D WITH_ONEBODY_OPT_AM="N0;N1;N2"`

* `--with-eri-max-am=N` --> `-D WITH_ERI_MAX_AM=N`
* `--with-eri-max-am=N0,N1,N2` --> `-D WITH_ERI_MAX_AM="N0;N1;N2"`
* `--with-eri-opt-am=N` --> `-D WITH_ERI_OPT_AM=N`
* `--with-eri-opt-am=N0,N1,N2` --> `-D WITH_ERI_OPT_AM="N0;N1;N2"`

* `--with-eri3-max-am=N` --> `-D WITH_ERI3_MAX_AM=N`
* `--with-eri3-max-am=N0,N1,N2` --> `-D WITH_ERI3_MAX_AM="N0;N1;N2"`
* `--with-eri3-opt-am=N` --> `-D WITH_ERI3_OPT_AM=N`
* `--with-eri3-opt-am=N0,N1,N2` --> `-D WITH_ERI3_OPT_AM="N0;N1;N2"`

* `--with-eri2-max-am=N` --> `-D WITH_ERI2_MAX_AM=N`
* `--with-eri2-max-am=N0,N1,N2` --> `-D WITH_ERI2_MAX_AM="N0;N1;N2"`
* `--with-eri2-opt-am=N` --> `-D WITH_ERI2_OPT_AM=N`
* `--with-eri2-opt-am=N0,N1,N2` --> `-D WITH_ERI2_OPT_AM="N0;N1;N2"`

* `--enable-shared` --> `-D BUILD_SHARED=ON` --> `-D BUILD_SHARED_LIBS=ON` (standard CMake variable)
* `--enable-static` --> `-D BUILD_STATIC=ON` --> `-D BUILD_SHARED_LIBS=OFF` (standard CMake variable)
* `--enable-shared --enable-static` --> `-D BUILD_SHARED=ON -D BUILD_STATIC=ON` --> `-D LIBINT2_BUILD_SHARED_AND_STATIC_LIBS=ON`
* `-D ENABLE_CXX11API=ON` --> `-D REQUIRE_CXX_API=ON`

* `--prefix=path` --> `-D CMAKE_INSTALL_PREFIX=path` (standard CMake variable)

* Targets
  * `libint2` --> `Libint2::int2` (internal target name `int-{shared,static}`)
  * `libint2_cxx` --> `Libint2::cxx` (internal target name `int-cxx-headeronly-{shared,static}`)
  * DNE --> `Libint2::int2-cxx` (internal target name `int-cxx-{shared,static}`)


