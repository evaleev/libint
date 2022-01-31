# libint compiler vs library

Before you read on:

* If you want to know how to _use_ a libint library in your code:
  * if you use C++11 or later (strongly recommended): read [this](https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API) instead
  * if you use pre-2011 C++, C, Fortran, or any other language, refer to the [Libint Programmer's Manual](http://sourceforge.net/projects/libint/files/libint-for-beginners/progman.pdf/download)
* If you want to know how to _generate_ a libint _library_ using the libint _compiler_, first make sure you really need to do that:
  * if all you want is a basic library that computes integrals necessary to compute energies, use the pre-generated library labeled "lmax=6 library (standard ints only)" from the [latest release](https://github.com/evaleev/libint/releases/latest) of Libint
  * many codes using libint, e.g. orca and mpqc, already include an appropriately configured libint library, and you do not need to generate it yourself
  * if you _do_ need to make a custom library, _read on_



# Prerequisites

| Task                                                                 | Compilers               | CMake[^3] | CMake generator | Boost[^7] | Eigen   | GMPXX[^13] | MPFR[^14] |
| :------------------------------------------------------------------- | :---------------------: | :-------: | --------------- | :-------: | :-----: | :--------: | :-------: |
| build target `build_libint`                                          | C++[^1]                 | 🔵[^4]    | Ninja           | 🔵[^8]    | &ndash; | 🔵         | &ndash;   |
| build target `library`                                               | C++[^1]                 | 🔵[^5]    | Ninja           | &ndash;   | &ndash; | &ndash;    | &ndash;   |
| &emsp;&emsp;`-D REQUIRE_CXX_API=ON`                                  | C++[^1]                 | 🔵[^5]    | Ninja           | 🔸[^9]    | 🔵[^11] | &ndash;    | &ndash;   |
| &emsp;&emsp;`-D ENABLE_FORTRAN=ON`                                   | C++[^1], Fortran[^2], C | 🔵[^5]    | Ninja           | &ndash;   | &ndash; | &ndash;    | &ndash;   |
| build&nbsp;project&nbsp;_consuming_&nbsp;Libint2&nbsp;library        |
| &emsp;C&nbsp;interface&nbsp;(I/F),&nbsp;`Libint2::int2`              | C++[^1]                 | 🔸[^6]    | Ninja, Makefile | &ndash;   | &ndash; | &ndash;    | &ndash;   |
| &emsp;C++11&nbsp;header&nbsp;I/F,&nbsp;`Libint2::cxx`                | C++[^1]                 | 🔸[^6]    | Ninja, Makefile | 🔸[^10]   | 🔵      | &ndash;    | &ndash;   |
| &emsp;&emsp;`-D ENABLE_MPFR=ON`                                      | C++[^1]                 | 🔸[^6]    | Ninja, Makefile | 🔸[^10]   | 🔵      | 🔵         | 🔵        |
| &emsp;C++11&nbsp;compiled&nbsp;I/F,&nbsp;`int2-cxx`                  | C++[^1]                 | 🔸[^6]    | Ninja, Makefile | 🔸[^10]   | 🔵[^12] | &ndash;    | &ndash;   |
| &emsp;Fortran I/F,&nbsp;`Libint2::fortran`                           | Fortran[^2]             | 🔸[^6]    | Ninja, Makefile |           |         | &ndash;    | &ndash;   |

* `🔵` required
* `🔸` required or recommended, but there's a path forward without
* `—` not involved

[^1]: C++ compiler that supports C++11 standard. C++11 standard is the fourth most recent international standard for C++, hence most modern compilers support it fully. A common compiler flag is `-std=c++11`, which CMake will impose on the compilation.

[^2]: Fortran 2003 compiler.

[^3]: CMake 3.16 or higher.

[^4]: Since Libint2 v2.8 TODO, the GNU toolchain has been replaced by CMake as the sole buildsystem for the Libint2 compiler, `build_libint`. See [update guide](#GNU-Autotools-Update-Guide).

[^5]: Since Libint2 v2.8 TODO, the CMake buildsystem for the exported library has been reworked. See [update guide](#GNU-Autotools-Update-Guide).

[^6]: Consuming an installed Libint2 library is simplest with CMake by employing `find_package(Libint2)` and `target_link_libraries(... Libint2::...)` commands. To facilitate consumption outside CMake, pkgconfig files are available for the C interface, and more could be provided.

[^7]: Boost 1.57 or higher. Only header-only (no compiled libraries) components needed.

[^8]: Building the Libint2 compiler needs several Boost components including MPL, Type Traits, and Preprocessor. A detectable system installation is required.

[^9]: Building the Libint2 library with C++11 API needs the Boost Preprocessor (PP) component. For the compiled C++11 interface, `Libint2::int2-cxx`, the PP is actually compiled against, but for the header-only target, `Libint2::cxx`, the PP only sets up the usage dependency. A system installation of Boost is sought, but if none suitable found, a bundled version of PP is installed within the Libint2 header namespace.

[^10]: Consuming an installed Libint2 library through a C++11 interface requires the Boost Preprocessor (PP) component. Depending on the library *build* environment, a copy may have been bundled/vendored with the install at `CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_INCLUDEDIR/libint2/boost/`.

[^11]: Building the Libint2 library with C++11 API needs the header-only Eigen library. For the compiled C++11 interface, `Libint2::cxx`, Eigen is actually compiled against, but for the header-only target `Libint2::cxx_ho`, Eigen only sets up the usage dependency. A detectable (either through Eigen3Config.cmake or through location-hinting) system installation is required.

[^12]: Consuming an installed Libint2 library through the compiled C++11 interface, `Libint2::int2-cxx` requires Eigen. It is *strongly* recommended that the same installation of Eigen be used both to build and consume the `Libint2::int2-cxx` target, especially as regards configuring BLAS and other backends.

[^13]: Building the Libint2 compiler or building the Libint2 library with `-D ENABLE_MPFR=ON` for high-precision testing requires the [GNU Multiple Precision (GMP)](https://gmplib.org/) library. A detectable system installation is required, and it must include C++ support. For Windows, the [MPIR](https://www.mpir.org) project may satisfy the requirement.

[^14]: Building against the Libint2 library for the purpose of high-precision testing with define `LIBINT_HAS_MPFR=1` requires the [MPFR](https://www.mpfr.org/) library. A detectable system installation is required.

# Synopsis

- configure: `cmake -S /path/to/compiler/source/tree -B build [-Dvar1=value1] [-Dvar2=value2] ...` where
  the optional CMake cache variables `var1`, `var2` are described below. Replace `build` with desired path to
  the build directory, if desired; it will be created, if needed.
- build: `cmake --build build`
- test (optional): `cmake --build build --target check`
- install the library: `cmake --build build --target install`


For (more) complete instructions please refer to https://github.com/evaleev/libint/wiki




CMake build overview:

```bash
>>> ls
cmake/  COPYING  src/  tests/  ...
>>> cmake -S. -Bbuild -GNinja -DCMAKE_INSTALL_PREFIX=/path/to/future/install-libint ...
...
-- Generating done
-- Build files have been written to: /current/dir/build
>>> cmake --build build --target install -j`getconf _NPROCESSORS_ONLN`
```


| target                    | incl. | steps |     ( |  see  | below | )     |
| ------                    | ----- | ----- | ----- | ----- | ----- | ----- |
| build_libint              |   1   |   -   |   -   |   -   |   -   |   -   |
| check-libint2compiler     |   1   |   2   |   -   |   -   |   -   |   -   |
| export                    |   1   |   -   |   3   |   -   |   -   |   -   |
| library (default)         |   1   |   -   |   3   |   4   |   -   |   -   |
| check                     |   1   |   2   |   3   |   4   |   5   |   -   |
| install                   |   1   |   -   |   3   |   4   |   -   |   6   |


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
  - can be built as a subproject (FetchContent) or completely insulated (bare ExternalProject; default)
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
# configuring libint generator

These are the most useful configure options:

* Notes
  * Codes "G", "L", or "GL" for each option indicate whether it is consumed by the generator, the library, or both.
    * If your final target is the export tarball, use options that include the letter "G".
    * If you're building a library from an export tarball, use options that include the letter "L".
    * For a continuous generator->export->library build, options supplied at the top level will be properly handed off to generator and library build.

###  Which Integrals Classes, Which Derivative Levels

* `ENABLE_ONEBODY` — G — Compile with support for up to N-th derivatives of 1-body integrals. Use -1 for OFF. [Default=0]
* `ENABLE_ERI` — G — Compile with support for up to N-th derivatives of 4-center electron repulsion integrals. Use -1 for OFF. [Default=0]
* `ENABLE_ERI3` — G — Compile with support for up to N-th derivatives of 3-center electron repulsion integrals. Use -1 for OFF. [Default=-1]
* `ENABLE_ERI2` — G — Compile with support for up to N-th derivatives of 2-center electron repulsion integrals. Use -1 for OFF. [Default=-1]
* `ENABLE_G12` — G — Compile with support for N-th derivatives of MP2-F12 energies with Gaussian factors. Use -1 for OFF. [Default=-1]
* `ENABLE_G12DKH` — G — Compile with support for N-th derivatives of DKH-MP2-F12 energies with Gaussian factors. Use -1 for OFF. [Default=-1]

* `DISABLE_ONEBODY_PROPERTY_DERIVS` — G — Disable geometric derivatives of 1-body property integrals (all but overlap, kinetic, elecpot).
   These derivatives are disabled by default to save compile time. Use OFF to enable. [Default=ON]

###  Which Ordering Convention

* `LIBINT2_SHGAUSS_ORDERING` — L — Ordering for shells of solid harmonic Gaussians. Consumed at library build-time. [Default=standard]
  * `standard` — standard ordering (-l, -l+1 ... l)
  * `gaussian` — the Gaussian ordering (0, 1, -1, 2, -2, ... l, -l)
* `LIBINT2_CARTGAUSS_ORDERING` — G — Orderings for shells of cartesian Gaussians. Consumed at generator build-time. [Default=standard]
  * `standard` — standard ordering (xxx, xxy, xxz, xyy, xyz, xzz, yyy, ...) This is ordering of the Common Component Architecture (CCA) standard for molecular integral data exchange described in ["Components for Integral Evaluation in Quantum Chemistry", J. P. Kenny, C. L. Janssen, E. F. Valeev, and T. L. Windus, J. Comp. Chem. 29, 562 (2008)](http://dx.doi.org/10.1002/jcc.20815).
  * `intv3` — intv3 ordering (yyy, yyz, yzz, zzz, xyy, xyz, xzz, xxy, xxz, xxx) This is used by IntV3, the default integral engine of [MPQC](https://github.com/evaleev/libint/wiki/www.mpqc.org). Use this to make Libint and IntV3 engines in MPQC interoperable.
  * `gamess` — [GAMESS](http://www.msg.ameslab.gov/gamess/) ordering (xxx, yyy, zzz, xxy, xxz, yyx, yyz, zzx, zzy, xyz)
  * `orca` — [ORCA](http://cec.mpg.de/forum/) ordering (hydrid between GAMESS and standard)
  * `bagel` — [BAGEL](https://github.com/evaleev/libint/wiki/nubakery.org) axis-permuted version of intv3 (xxx, xxy, xyy, yyy, xxz, xyz, yyz, xzz, yzz, zzz)
* `LIBINT2_SHELL_SET` — G — Support computation of shell sets sets subject to these restrictions. Consumed at generator build-time. [Default=standard]
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
* `ERI3_PURE_SH` — G — Assume the 'unpaired' center of 3-center ERIs will be transformed to pure solid harmonics. [Default=OFF]
* `ERI2_PURE_SH` — G — Assume the 2-center ERIs will be transformed to pure solid harmonics. [Default=OFF]

###  How High Angular Momentum

* Notes
  * example for "semicolon-separated string": `-DENABLE_ERI3=2 -DWITH_ERI3_MAX_AM="5;4;3"`. cmake configuration prints:

    ```
    -- Setting option ENABLE_ERI3: 2
    -- Setting option WITH_ERI3_MAX_AM: 5;4;3
    ```

  * high MAX_AM generating >20k files may require `ulimit -s 65535` for linking library target on Linux to avert "ld: Argument list too long"

* `WITH_MAX_AM` — G — Support Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=4]
* `WITH_OPT_AM` — G — Optimize maximally for up to angular momentum N (N <= WITH_MAX_AM). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `(WITH_MAX_AM/2)+1`]

* `MULTIPOLE_MAX_ORDER` — G — Maximum order of spherical multipole integrals. There is no maximum. [Default=4]

* `WITH_ONEBODY_MAX_AM` — G — Support 1-body ints for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ONEBODY_OPT_AM` — G — Optimize 1-body ints maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI_MAX_AM` — G — Support 4-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI_OPT_AM` — G — Optimize 4-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI3_MAX_AM` — G — Support 3-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI3_OPT_AM` — G — Optimize 3-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI2_MAX_AM` — G — Support 2-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI2_OPT_AM` — G — Optimize 2-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

### Miscellaneous

* Notes
  * Approximate defaults are shown. Actual defaults from [GNUInstallDirs](https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html)

* `CMAKE_INSTALL_BINDIR` — L — Directory to which executables and runtime libraries are installed. Standard CMake variable. [Default=bin]
* `CMAKE_INSTALL_LIBDIR` — L — Directory to which libraries are installed. Standard CMake variable. [Default=lib]
* `CMAKE_INSTALL_INCLUDEDIR` — L — Directory to which headers are installed. Standard CMake variable. [Default=include]
* `CMAKE_INSTALL_DATADIR` — L — Directory to which data files are installed. Standard CMake variable. [Default=share]
* `LIBINT2_INSTALL_CMAKEDIR` — L — Directory to which CMake files are installed. [Default=lib/cmake/libint2]
* `LIBINT2_INSTALL_BASISDIR` — L — Directory to which data (basis) files are installed. basis/ directory created within this. [Default=share/libint/<LIBINT_VERSION>]

* `BUILD_TESTING` — GL — Whether to build the testing infrastructure and define the `check` target. Standard CMake variable. [Default=ON]
* `ENABLE_MPFR` — L — Use MPFR dependency to test Libint integrals in high precision. [Default=OFF]

## GNU Autotools Update Guide

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

* `--with-multipole-max-order=N` --> `-D MULTIPOLE_MAX_ORDER=N`

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

* `--enable-mpfr` --> assumed present --> `-D ENABLE_MPFR=ON`

* `--prefix=path` --> `-D CMAKE_INSTALL_PREFIX=path` (standard CMake variable)

* Targets
  * All targets listed below are available. rhs of arrow targets (namespaced with `Libint2::`) are preferred. lhs of arrow targets have legacy aliases for now.
  * Namespaced targets available through `find_package()` or `add_subdirectory()`
  * `libint2` --> `Libint2::int2` (internal target name `int-{shared,static}`)
  * `libint2_cxx` --> `Libint2::cxx` (internal target name `int-cxx-headeronly-{shared,static}`)
  * DNE --> `Libint2::int2-cxx` (internal target name `int-cxx-{shared,static}`)


## Packagers

* Decide if you want the Boost preprocessor headers bundled with Libint or if they should be a
  build-against-time dependency of the C++11 interface. Withhold (bundle) or supply (dependency)
  Boost detection paths from the library build accordingly. FWIW, Conda bundles.

## program-specific notes

#### mpqc4

* standard libtool configuration:

  ```
  --enable-generic-code --with-max-am=6 --with-opt-am=3 --enable-eri3=0 --enable-eri2=0 --enable-eri3-pure-sh --enable-eri2-pure-sh --enable-fma --disable-1body-property-derivs
  ```

* libtool configuration prior to Jan. 8, 2015:

  ```
  --enable-eri=0 --with-max-am=7 --with-opt-am=4 --disable-unrolling --enable-generic-code --enable-contracted-ints
  ```

#### gamess

* standard libtool configuration:

  ```
  --enable-eri=0 --with-max-am=7 --with-opt-am=4 --disable-unrolling --enable-generic-code --enable-contracted-ints --with-cartgauss-ordering=gamess
  ```

#### orca

* a libint library (version 2.0.2) is embedded in ORCA
* standard libtool configuration:

  ```
  --enable-eri=2 --enable-eri3=2 --enable-eri2=2 --with-max-am=7 --with-opt-am=4 --with-eri-max-am=7,4,3 --with-eri-opt-am=4,3,2 --disable-unrolling --enable-generic-code --enable-contracted-ints --with-cartgauss-ordering=orca --with-shell-set=orca --enable-eri3-pure-sh --enable-eri2-pure-sh
  ```

#### bagel

* standard libtool configuration:

  ```
  --with-max-am=4 --with-eri3-max-am=6 --with-eri2-max-am=6 --enable-eri3=1 -enable-eri=1 --enable-eri2=1 --disable-unrolling --enable-generic-code --enable-contracted-ints --with-cartgauss-ordering=bagel

  ```

* if you want to use spherical Gaussians only add: `--enable-eri3-pure-sh --enable-eri2-pure-sh` (some tests may fail)
* It appears that on a Mac Libint and BAGEL must be either both static or both shared (2/3/2014)

#### psi4

* production CMake configuration:

  ```
  -DREQUIRE_CXX_API=ON -DLIBINT2_SHGAUSS_ORDERING=gaussian -DLIBINT2_CARTGAUSS_ORDERING=standard -DLIBINT2_SHELL_SET=standard -DENABLE_ERI=2 -DENABLE_ERI3=2 -DENABLE_ERI2=2 -DENABLE_ONEBODY=2 -DMULTIPOLE_MAX_ORDER=4 -DWITH_MAX_AM="6;5;4" -DWITH_ERI_MAX_AM="5;4;3" -DWITH_ERI3_MAX_AM="6;5;4" -DWITH_ERI2_MAX_AM="6;5;4" -DERI3_PURE_SH=OFF -DERI2_PURE_SH=OFF
  ```

* minimal detection:

  ```
  find_package(Libint2 COMPONENTS gss CXX_ho impure_sh eri_c4_d0_l3 eri_c4_d1_l2 eri_c4_d2_l2 eri_c3_d0_l4 eri_c3_d1_l3 eri_c3_d2_l3 eri_c2_d0_l4 eri_c2_d1_l3 eri_c2_d2_l3 onebody_d0_l4 onebody_d1_l3 onebody_d2_l3)
  ```

* see [notes](https://github.com/psi4/psi4/blob/master/external/upstream/libint2/CMakeLists.txt) for details of Psi4/Libint2 configuration

    -DDISABLE_ONEBODY_PROPERTY_DERIVS=OFF \
    -DENABLE_G12=1 \
    -DWITH_G12_MAX_AM=4 \
    -DWITH_G12_OPT_AM=3 \

