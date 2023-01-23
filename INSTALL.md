# Libint Compiler vs Library

Before you read on:

* If you want a pre-built libint library, packages may be available:
  * conda-forge (TBD)
  * Debian (TBD)
  * Fedora (TBD)
* If you want to know how to _use_ a libint library in your code:
  * if you use C++11 or later (strongly recommended): read [this](https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API) instead
  * if you use pre-2011 C++, C, Fortran, or any other language, refer to the [Libint Programmer's Manual](http://sourceforge.net/projects/libint/files/libint-for-beginners/progman.pdf/download)
* If you want to _build and use_ a libint _library_:
  * if all you want is a basic library that computes integrals necessary to compute energies, use the pre-generated library labeled "lmax=6 library (standard ints only)" from the [latest release](https://github.com/evaleev/libint/releases/latest) of Libint
  * many codes using libint, e.g. orca and mpqc, already include an appropriately configured libint library, and you do not need to generate it yourself
  * if you need compilation directions, _read on_, skipping the compiler/generation parts.
* If you want to know how to _generate_ a libint _library_ using the libint _compiler_, these are some compelling circumstances:
  * if you need a custom libint library with choice of integral types, AM, orderings, language interfaces, etc.
  * if you want to develop libint with new integral types, recurrence relations, and computation
    strategies, you'll need to edit the compiler. If you are interested in working on the compiler
    code please consider consulting with one of the Libint authors to avoid duplication of effort.
  * if you do need to generate a custom library, _read on_.


-----------------------------------------------------------------------------

# Overview

The Libint build is structured into three parts:

* generator/compiler
  - (1) build src/bin/libint/ into compiler executable `build_libint`
  - pretty quick, runs in parallel
  - consumes all the enable/max/opt integral options and all orderings options except solid harmonic
  - (2) optionally testable
* export
  - (3) run `build_libint` to generate library source (C++) files (that upon
    compilation can become a Libint2 library) and combine them with other
    static source files in src/lib/libint/ and general static files (e.g.,
    include/ and docs/) into an independent tarball ready for distribution
    (with its own CMake configuration, tests, etc.).
  - really slow for non-trivial angular momenta, runs in serial
  - consumes no options
  - build target `export` to stop after this step and collect source tarball
* library
  - can be built as a subproject (FetchContent) or completely insulated (bare ExternalProject; default).
    For FetchContent, must build libint-library-export target before library build targets appear
  - (4) unpack the export tarball and build the library and install into \<build\>/library-install-stage/
  - duration depends on number of integrals requested, runs in parallel
  - consumes solid harmonic ordering and the CMAKE_INSTALL_[DATA|INCLUDE|LIB]DIR
  - the default build target includes this final library build
  - (5) optionally testable
  - (6) install into CMAKE_INSTALL_PREFIX


Command-line synopsis. See [table](#Build-Targets) for `--target` choices (steps refer to numbered bullets above) and [section](#Configuring-Libint) for `-D options` choices.

```bash
>>> git clone https://github.com/evaleev/libint.git && cd libint
>>> ls
cmake/  COPYING  src/  tests/  ...
>>> cmake -S. -Bbuild -GNinja -DCMAKE_INSTALL_PREFIX=/path/to/future/install-libint -D options ...
...
-- Generating done
-- Build files have been written to: /current/dir/build
>>> cmake --build build --target install -j`getconf _NPROCESSORS_ONLN`
```

### Build Targets

| `--target ...`            | incl. | steps |     ( |  see  | above | )     |
| --------------            | ----- | ----- | ----- | ----- | ----- | ----- |
| `build_libint`            |   1   |   -   |   -   |   -   |   -   |   -   |
| `check-libint2compiler`   |   1   |   2   |   -   |   -   |   -   |   -   |
| `export`                  |   1   |   -   |   3   |   -   |   -   |   -   |
| `library` (default)       |   1   |   -   |   3   |   4   |   -   |   -   |
| `check`                   |   1   |   2   |   3   |   4   |   5   |   -   |
| `install`                 |   1   |   -   |   3   |   4   |   -   |   6   |
| `check install`           |   1   |   2   |   3   |   4   |   5   |   6   |

Use combined targets like `cmake --target check install` to avoid some unnecessary rebuilding (esp. of build_libint) that occurs with successive targets. The CMake dependency structure is imperfect.


-----------------------------------------------------------------------------

# Prerequisites

| Task                                                                 | Compilers               | CMake[^3] | CMake generator[^20] | Py      | Boost[^7] | Eigen   | GMPXX[^13] | MPFR[^14] |
| :------------------------------------------------------------------- | :---------------------: | :-------: | -------------------- | :-----: | :-------: | :-----: | :--------: | :-------: |
| build target `build_libint`                                          | C++[^1]                 | 🔵[^4]    | Ninja                | &ndash; | 🔵[^8]    | &ndash; | 🔵         | &ndash;   |
| build target `library`                                               | C++[^1], C              | 🔵[^5]    | Ninja                | 🔸[^21] | &ndash;   | &ndash; | &ndash;    | &ndash;   |
| &emsp;&emsp;`-D REQUIRE_CXX_API=ON`                                  | C++[^1], C              | 🔵[^5]    | Ninja                | 🔸[^21] | 🔸[^9]    | 🔵[^11] | &ndash;    | &ndash;   |
| &emsp;&emsp;`-D ENABLE_FORTRAN=ON`                                   | C++[^1], Fortran[^2], C | 🔵[^5]    | Ninja                | 🔵[^22] | &ndash;   | &ndash; | &ndash;    | &ndash;   |
| &emsp;&emsp;`-D ENABLE_PYTHON=ON`                                    | C++[^1], C              | 🔵[^5]    | Ninja                | 🔵[^23] | 🔸[^9]    | 🔵[^11] | &ndash;    | &ndash;   |
| build&nbsp;project&nbsp;_consuming_&nbsp;Libint2&nbsp;library        |
| &emsp;C&nbsp;interface&nbsp;(I/F),&nbsp;`Libint2::int2`              | C++[^1]                 | 🔸[^6]    | Ninja, Makefile      | &ndash; | &ndash;   | &ndash; | &ndash;    | &ndash;   |
| &emsp;C++11&nbsp;header&nbsp;I/F,&nbsp;`Libint2::cxx`                | C++[^1]                 | 🔸[^6]    | Ninja, Makefile      | &ndash; | 🔸[^10]   | 🔵      | &ndash;    | &ndash;   |
| &emsp;&emsp;`-D ENABLE_MPFR=ON`                                      | C++[^1]                 | 🔸[^6]    | Ninja, Makefile      | &ndash; | 🔸[^10]   | 🔵      | 🔵         | 🔵        |
| &emsp;C++11&nbsp;compiled&nbsp;I/F,&nbsp;`int2-cxx`                  | C++[^1]                 | 🔸[^6]    | Ninja, Makefile      | &ndash; | 🔸[^10]   | 🔵[^12] | &ndash;    | &ndash;   |
| &emsp;Fortran I/F,&nbsp;`Libint2::fortran`                           | Fortran[^2]             | 🔸[^6]    | Ninja, Makefile      | &ndash; |           |         | &ndash;    | &ndash;   |

* `🔵` required
* `🔸` required or recommended, but there's a path forward without
* `—` not involved

[^1]: C++ compiler that supports C++11 standard. C++11 standard is the fourth most recent international standard for C++, hence most modern compilers support it fully. A common compiler flag is `-std=c++11`, which CMake will impose on the compilation.

[^2]: Fortran 2003 compiler to enable Fortran bindings generation.

[^3]: [CMake](https://cmake.org/) 3.16 or higher.

[^4]: Since Libint2 v2.8 TODO, the GNU toolchain has been replaced by CMake as the sole buildsystem for the Libint2 compiler, `build_libint`. See [update guide](#GNU-Autotools-Update-Guide).

[^5]: Since Libint2 v2.8 TODO, the CMake buildsystem for the exported library has been reworked. See [update guide](#GNU-Autotools-Update-Guide).

[^6]: Consuming an installed Libint2 library is simplest with CMake by employing `find_package(Libint2)` and `target_link_libraries(... Libint2::...)` commands. To facilitate consumption outside CMake, pkgconfig files are available for the C interface, and more could be provided.

[^7]: [Boost](https://www.boost.org/) 1.57 or higher. Only header-only (no compiled libraries) components needed.

[^8]: Building the Libint2 compiler needs several Boost components including MPL, Type Traits, and Preprocessor. A detectable system installation is required. (That is, "bundled Boost" is insufficient.)

[^9]: Building the Libint2 library with C++11 API needs the Boost Preprocessor (PP) component. For the compiled C++11 interface, `Libint2::int2-cxx`, the PP is actually compiled against, but for the header-only target, `Libint2::cxx`, the PP only sets up the usage dependency. A system installation of Boost is sought, but if none suitable found, a bundled version of PP is installed within the Libint2 header namespace.

[^10]: Consuming an installed Libint2 library through a C++11 interface requires the Boost Preprocessor (PP) component. Depending on the library *build* environment, a copy may have been bundled/vendored with the install at `CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_INCLUDEDIR/libint2/boost/`.

[^11]: Building the Libint2 library with C++11 API needs the header-only [Eigen](https://eigen.tuxfamily.org/) library. For the compiled C++11 interface, `Libint2::cxx`, Eigen is actually compiled against, but for the header-only target `Libint2::cxx_ho`, Eigen only sets up the usage dependency. A detectable (either through Eigen3Config.cmake or through location-hinting) system installation is required.

[^12]: Consuming an installed Libint2 library through the compiled C++11 interface, `Libint2::int2-cxx` requires [Eigen](https://eigen.tuxfamily.org/). It is *strongly* recommended that the same installation of Eigen be used both to build and consume the `Libint2::int2-cxx` target, especially as regards configuring BLAS and other backends.

[^13]: Building the Libint2 compiler or building the Libint2 library with `-D ENABLE_MPFR=ON` for high-precision testing requires the [GNU Multiple Precision (GMP)](https://gmplib.org/) library. A detectable system installation is required, and it must include C++ support. For Windows, the [MPIR](https://www.mpir.org) project satisfies the requirement.

[^14]: Building against the Libint2 library for the purpose of high-precision testing with define `LIBINT_HAS_MPFR=1` requires the [MPFR](https://www.mpfr.org/) library. A detectable system installation is required.

[^20]: Tested CMake generators are [Ninja](https://ninja-build.org/) or [GNU Make](https://www.gnu.org/software/make/). The use of Ninja is **strongly** recommended!

[^21]: Python used for testing.

[^22]: Python used to process files for Fortran binding.

[^23]: Python headers and interpreter needed for Pybind11 module

-----------------------------------------------------------------------------

# Configuring Libint

* Notes
  * Codes "G", "L", or "C" for each option indicate whether it is consumed by the _g_enerator, the _l_ibrary, the library _c_onsumer, or a combination.
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

* `ENABLE_T1G12_SUPPORT` — G — Enable [Ti,G12] integrals when G12 integrals are enabled. Irrelevant when `ENABLE_G12=OFF`. Use OFF to disable. [Default=ON]


###  Which Ordering Conventions

* `LIBINT2_SHGAUSS_ORDERING` — L — Ordering for shells of solid harmonic Gaussians. [Default=standard]
  * `standard` — standard ordering (-l, -l+1 ... l)
  * `gaussian` — the Gaussian ordering (0, 1, -1, 2, -2, ... l, -l)
* `LIBINT2_CARTGAUSS_ORDERING` — G — Orderings for shells of cartesian Gaussians. [Default=standard]
  * `standard` — standard ordering (xxx, xxy, xxz, xyy, xyz, xzz, yyy, ...) This is ordering of the Common Component Architecture (CCA) standard for molecular integral data exchange described in ["Components for Integral Evaluation in Quantum Chemistry", J. P. Kenny, C. L. Janssen, E. F. Valeev, and T. L. Windus, J. Comp. Chem. 29, 562 (2008)](http://dx.doi.org/10.1002/jcc.20815).
  * `intv3` — intv3 ordering (yyy, yyz, yzz, zzz, xyy, xyz, xzz, xxy, xxz, xxx) This is used by IntV3, the default integral engine of [MPQC](https://github.com/evaleev/libint/wiki/www.mpqc.org). Use this to make Libint and IntV3 engines in MPQC interoperable.
  * `gamess` — [GAMESS](http://www.msg.ameslab.gov/gamess/) ordering (xxx, yyy, zzz, xxy, xxz, yyx, yyz, zzx, zzy, xyz)
  * `orca` — [ORCA](http://cec.mpg.de/forum/) ordering (hydrid between GAMESS and standard)
  * `bagel` — [BAGEL](https://github.com/evaleev/libint/wiki/nubakery.org) axis-permuted version of intv3 (xxx, xxy, xyy, yyy, xxz, xyz, yyz, xzz, yzz, zzz)
* `LIBINT2_SHELL_SET` — G — Support computation of shell sets sets subject to these restrictions. [Default=standard]
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

  * high MAX_AM generating >20k files may require `ulimit -s 65535` for linking library target on Linux to avert "ld: Argument list too long". Unity build (ON by default for library) averts this.

* `WITH_MAX_AM` — G — Support Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. Specify values greater or equal to `WITH_<class>_MAX_AM`; often mirrors `WITH_ERI3_MAX_AM`. [Default=4]
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

* `WITH_G12_MAX_AM` — G — Support integrals for G12 methods of angular momentum up to N. No specification with per-derivative list. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_G12_OPT_AM` — G — Optimize G12 integrals for up to angular momentum N (N <= max-am). No specification with per-derivative list. [Default=-1 `WITH_OPT_AM`]

* `WITH_G12DKH_MAX_AM` — G — Support integrals for relativistic G12 methods of angular momentum up to N. No specification with per-derivative list. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_G12DKH_OPT_AM` — G — Optimize G12DKH integrals for up to angular momentum N (N <= max-am). No specification with per-derivative list. [Default=-1 `WITH_OPT_AM`]


### Compilers and Flags

* `CMAKE_CXX_COMPILER` — G L — Specify C++ Compiler. [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER.html)
* `CMAKE_C_COMPILER` — G L — Specify C Compiler. [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER.html)
* `CMAKE_Fortran_COMPILER` — G L — Specify Fortran Compiler. [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER.html)
* `CMAKE_CXX_FLAGS` — G L — Additional C++ flags. [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_FLAGS.html)
* `CMAKE_C_FLAGS` — G L — Additional C flags. [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_FLAGS.html)
* `CMAKE_Fortran_FLAGS` — G L — Additional Fortran flags. [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_FLAGS.html)


### Install Paths

* Notes
  * Approximate defaults are shown. Actual defaults from [GNUInstallDirs](https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html)

* `CMAKE_INSTALL_PREFIX` — L — Directory into which library installed. [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html)
* `CMAKE_INSTALL_BINDIR` — L — Directory within `CMAKE_INSTALL_PREFIX` to which executables and runtime libraries are installed. [Standard CMake variable](https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html) [Default=bin]
* `CMAKE_INSTALL_LIBDIR` — L — Directory within `CMAKE_INSTALL_PREFIX` to which libraries are installed. [Standard CMake variable](https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html) [Default=lib]
* `CMAKE_INSTALL_INCLUDEDIR` — L — Directory within `CMAKE_INSTALL_PREFIX` to which headers are installed. [Standard CMake variable](https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html) [Default=include]
* `CMAKE_INSTALL_DATADIR` — L — Directory within `CMAKE_INSTALL_PREFIX` to which data files are installed. [Standard CMake variable](https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html) [Default=share]
* `LIBINT2_INSTALL_CMAKEDIR` — L — Directory within `CMAKE_INSTALL_PREFIX` to which CMake files are installed. [Default=lib/cmake/libint2]
* `LIBINT2_INSTALL_BASISDIR` — L — Directory within `CMAKE_INSTALL_PREFIX` to which data (basis) files are installed. basis/ directory created within this. [Default=share/libint/<LIBINT_VERSION>]
* `LIBINT2_INSTALL_FMODDIR` — L — Directory within `CMAKE_INSTALL_PREFIX` to which Fortran module files are installed if `ENABLE_FORTRAN=ON`. [Default=include/libint2/fortran2/modules]
* `PREFIX_PYTHON_INSTALL` — L — For `ENABLE_PYTHON=ON`, whether to install the Python module in the Linux manner to `CMAKE_INSTALL_PREFIX` or to not install it. Note: not a path; the installation sub-path below `CMAKE_INSTALL_PREFIX` is determined by querying `Python_EXECUTABLE`. For alternate installation in the Python manner to `Python_EXECUTABLE`'s site-packages, see target libint2-python-wheel. [Default=OFF]


### Detecting Dependencies

* `Python_EXECUTABLE` — L — Path to Python interpreter.
* `CMAKE_PREFIX_PATH` — G L — Set to list of root directories to look for external dependencies. [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/CMAKE_PREFIX_PATH.html)
* `BOOST_ROOT` — G L C —
* `Multiprecision_ROOT` — G L —
* `Eigen3_ROOT` — L C — Prefix to installation location (`Eigen3_ROOT/...` exists)
* `Libint2_DIR` — C — CMake variable, set to directory containing this Config file
* `LIBINT_LOCAL_Eigen3_FIND` — C — Set to `ON` before `find_package(Libint2)` to load the Eigen3 target exported by `LIBINT_LOCAL_Eigen3_INSTALL=ON` if Libint library built locally. [Default=OFF]
* `CMAKE_DISABLE_FIND_PACKAGE_Boost` — L — When Boost required for C++11 Libint API, disable its detection, thereby forcing use of bundled Boost. Note that this (and other Boost-hinting variables) can affect what is installed [see here](#packagers). [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/CMAKE_DISABLE_FIND_PACKAGE_PackageName.html). [Default=OFF]

EIGEN3_INCLUDE_DIR?

* Hint dependency locations all at the same installation prefix:

  ```
  -D CMAKE_PREFIX_PATH="/path/to/installation/prefix"
  -D CMAKE_PREFIX_PATH="/home/miniconda/envs/l2dev"
  ```

* Hint dependency locations all at different installation prefixes:

  ```
  -D CMAKE_PREFIX_PATH="/home/miniconda/envs/onlyboost;/home/miniconda/envs/onlygmp;/home/miniconda/envs/onlyeigen"
  ```

* Hint dependency locations targeted by package:

  ```
  -D BOOST_ROOT="/home/miniconda/envs/onlyboost"
  -D Multiprecision_ROOT="/home/miniconda/envs/onlygmp"
  -D Eigen3_ROOT="/home/miniconda/envs/onlyeigen"
  ```

* Hint dependency locations targeted by <package>Config.cmake (most CMake-like):

  ```
  -D Eigen3_DIR="/home/miniconda/envs/onlyeigen/share/eigen3/cmake"
  -D Boost_DIR="/home/miniconda/envs/onlyboost/lib/cmake/Boost-1.73.0"
  ```

* Hint dependency locations targeted by package variables (least CMake-like): UNTESTED

  ```
  ```


### Build Library What

* `REQUIRE_CXX_API` — L — Build C++11 Libint API. Define header-only library target and check target (requires Eigen3; Boost recommended; [see prereq line](#prerequisites)). [Default=ON]
* `REQUIRE_CXX_API_COMPILED` — L — Build C++11 Libint API. Define compiled (not just header-only) targets (requires Eigen3; Boost recommended). [Default=ON]
* `ENABLE_FORTRAN` — L — Build Fortran03+ module/bindings (requires C and Fortran compilers and Python). [Default=OFF]
* `ENABLE_MPFR` — L — Use MPFR library to test Libint integrals in high precision (requires MPFR; experts only). [Default=OFF]
* `LIBINT_LOCAL_Eigen3_INSTALL` — L — Install an exported target with hard-coded Eigen3 dependency paths. This is potentially useful and important when consuming the compiled C++11 interface library so that the Libint library build and Libint consumer build use the same Eigen3 installation & ABI. This is at most a convenience when consuming the header-only C++11 interface library. See `LIBINT_LOCAL_Eigen3_FIND`. [Default=OFF]
* `ENABLE_PYTHON` — L — Build Python bindings (requires Python and Eigen3; Boost and pybind11 recommended; [see prereq line](#prerequisites)). Can instead be enabled and built through separate CMake configuration after library build. [Default=OFF]


### Build Library How

* `CMAKE_BUILD_TYPE` — G L — [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html) [Default=Release]
* `BUILD_SHARED_LIBS` — L — Build Libint library as shared, not static. [Standard CMake variable](https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html) [Default=OFF]
* `LIBINT2_BUILD_SHARED_AND_STATIC_LIBS` — L — Build both shared and static Libint libraries in one shot. Uses `-fPIC`. [Default=OFF]
* `ENABLE_XHOST` — L — Enables processor-specific optimization (with MSVC, it enables AVX2 instructions) [Default=ON]
* `BUILD_TESTING` — G L — Whether to build the testing infrastructure and define the `check` target. [Standard CMake variable](https://cmake.org/cmake/help/latest/command/enable_testing.html) [Default=ON]
* `LIBINT_BUILD_LIBRARY_AS_SUBPROJECT` — G — If building compiler and library in continuous command, build generated library as a subproject; if OFF will configure and build separately (expert only). [Default=OFF]


### Miscellaneous

* `LIBINT2_REALTYPE` — L — Specifies the floating-point data type used by the library. [Default=double]
  By overriding the default it is possible to customize the library to use a lower-precision representation (which typically results in a performance boost) and/or to generate [SIMD](http://en.wikipedia.org/wiki/SIMD) vectorized code. *N.B. C++11 interface cannot be currently used with SIMD vectorized libraries!* The following values are valid:
  * `double` -- double-precision floating-point representation of a real number;
  * `float` -- single-precision floating-point number;
  * `libint2::simd::VectorAVXDouble` -- vector of 4 packed doubles that can be used with [AVX](http://en.wikipedia.org/wiki/Advanced_Vector_Extensions) instructions available on reasonably-modern x86 hardware (starting with Intel Sandy Bridge and AMD Bulldozer microarchitectures, available in processors since 2011);
  * `libint2::simd::VectorSSEDouble` -- vector of 2 packed doubles that can be used with [SSE2](http://en.wikipedia.org/wiki/SSE2) instructions available on all x86 platforms, including those released before 2011;
  * `libint2::simd::VectorSSEFloat`  -- vector of 4 packed floats that can be used with [SSE](http://en.wikipedia.org/wiki/Streaming_SIMD_Extensions) instructions available on all x86 platforms, including those released before 2011;
  * `libint2::simd::VectorQPXDouble` -- vector of 4 packed doubles that can be used with QPX instructions available on recent PowerPC hardware (IBM Blue Gene/Q);
  * `libint2::simd::VectorFP2Double` -- vector of 2 packed doubles that can be used with FP2 (Double Hummer) instructions available on older PowerPC hardware (IBM Blue Gene/P).

  With the exception of `float`, these are vector types implemented in Libint using compiler _intrinsics_, functions that translate directly into vector instructions. To use these vector types you may need to provide additional compiler flags that will enable support for vector instructions. For example, to enable support for AVX in Clang use the `-mavx` compiler flag. With Intel compiler use flag `-xHOST` to enable all vector instruction sets supported by the processor on which you are compiling.

  **N.B.** It is also possible to use real vector types of [Agner Fog's vectorclass library](http://www.agner.org/optimize/#vectorclass), e.g. `Vec4d` and `Vec8f` for AVX. To use this library you need to add this to CPPFLAGS or CXXFLAGS: `-Ipath_to_vectorclass -DLIBINT2_HAVE_AGNER_VECTORCLASS` . On macOS, we only succeeded in using this library with a recent GNU C++ compiler, not with Clang. Not tested after CMake rework.

* `LIBINT_CONTRACTED_INTS` — G — Turn on support for contracted integrals. [Default=ON]
* `LIBINT_ERI_STRATEGY` — G — Compute ERIs using the following strategy (experts only). [Default=1]
* `LIBINT_USE_COMPOSITE_EVALUATORS` — G — Libint will use composite evaluators (i.e. every evaluator will compute one integral type only). [Default=ON]
* `LIBINT_SINGLE_EVALTYPE` — G — Generate single evaluator type (i.e. all tasks use the same evaluator). [Default=ON]
* `LIBINT_ENABLE_UNROLLING` — G — Unroll shell sets into integrals (will unroll shell sets larger than N) (no->0, yes->1000000000). [Default=100]
* `LIBINT_ALIGN_SIZE` — G — If posix_memalign is available, this will specify alignment of Libint data, in units of sizeof(LIBINT2_REALTYPE). Default is to use built-in heuristics (experts only). [Default=0]
* `LIBINT_GENERATE_FMA` — G — Generate FMA (fused multiply-add) instructions (to benefit must have FMA-capable hardware and compiler). [Default=OFF]
* `LIBINT_ENABLE_GENERIC_CODE` — G — Use manually-written generic code. [Default=OFF]


-----------------------------------------------------------------------------

# GNU Autotools Update Guide

* Notes
  * Multiple option names can be from any long-lived branch but usually libtool+cmake --> final cmake+cmake.

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

* `--enable-g12=N` --> `-D ENABLE_G12=N`
* `--enable-g12dkh=N` --> `-D ENABLE_G12DKH`
* `--disable-t1g12-support` --> `-D ENABLE_T1G12_SUPPORT=OFF`
* `--with-g12-max-am=N` --> `-D WITH_G12_MAX_AM=N`
* `--with-g12-opt-am=N` --> `-D WITH_G12_OPT_AM=N`
* `--with-g12dkh-max-am=N` --> `-D WITH_G12DKH_MAX_AM=N`
* `--with-g12dkh-opt-am=N` --> `-D WITH_G12DKH_OPT_AM=N`

* `--disable-1body-property-derivs` --> `-D DISABLE_ONEBODY_PROPERTY_DERIVS=ON`

* `--enable-shared` --> `-D BUILD_SHARED=ON` --> `-D BUILD_SHARED_LIBS=ON` (standard CMake variable)
* `--enable-static` --> `-D BUILD_STATIC=ON` --> `-D BUILD_SHARED_LIBS=OFF` (standard CMake variable)
* `--enable-shared --enable-static` --> `-D BUILD_SHARED=ON -D BUILD_STATIC=ON` --> `-D LIBINT2_BUILD_SHARED_AND_STATIC_LIBS=ON`

* `-D REQUIRE_CXX_API=ON` --> `-D ENABLE_CXX11API=ON` --> `-D REQUIRE_CXX_API=ON`
* `--enable-mpfr` --> assumed present --> `-D ENABLE_MPFR=ON`
* `--prefix=path` --> `-D CMAKE_INSTALL_PREFIX=path` (standard CMake variable)
* `--with-cmakedir=partialpath` --> `-D LIBINT2_INSTALL_CMAKEDIR=partialpath`
* `--with-real-type=type` --> `-D LIBINT2_REALTYPE=type`

* (target) `libint2` --> `Libint2::int2`
* (target) `libint2_cxx` --> `Libint2::cxx`

* `ENV(CXX)=/path/to/c++/compiler` --> `-D CMAKE_CXX_COMPILER=/path/to/c++/compiler`
* `ENV(CXXFLAGS)` --> `-D CMAKE_CXX_FLAGS`
* `ENV(CPPFLAGS)=-I/path/to/boost/includes` --> `-D BOOST_ROOT=/path/to/boost/prefix`
* `ENV(FC)=/path/to/fortran/compiler` --> `-D CMAKE_Fortran_COMPILER=/path/to/fortran/compiler`

* `-D LIBINT2_PYTHON=ON` --> `-D ENABLE_PYTHON=ON`
* `-D LIBINT_USE_BUNDLED_BOOST=ON` --> `-D CMAKE_DISABLE_FIND_PACKAGE_Boost=ON` (standard CMake variable)

-----------------------------------------------------------------------------

# Consuming Libint

### Programming to Access Integrals

* if you use C++11 or later (strongly recommended): read [this](https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API) instead
* if you use pre-2011 C++, C, Fortran, or any other language, refer to the [Libint Programmer's Manual](http://sourceforge.net/projects/libint/files/libint-for-beginners/progman.pdf/download) for brief information on how to use the library in your code.

### Consumption Targets

| Namespaced Target[^15] | Component[^16] | Built by Default | Ensure Built                  | Ensure Excluded                           | Internal Target(s)[^17]              | Alias[^18]    |
| ---------------------- | -------------- | ---------------- | ----------------------------- | ----------------------------------------- | -----------------------------------  | ------------  |
| `Libint2::int2`        | `C`            | yes              | always                        | impossible                                | `int-{static,shared}`                | `libint2`     |
| `Libint2::cxx`         | `CXX_ho`       | yes              | `REQUIRE_CXX_API=ON`          | `REQUIRE_CXX_API=OFF` and withhold Eigen3 | `int-cxx-headeronly-{static,shared}` | `libint2_cxx` |
| `Libint2::int2-cxx`    | `CXX`          | yes              | `REQUIRE_CXX_API_COMPILED=ON` | `REQUIRE_CXX_API_COMPILED=OFF`            | `int-cxx-compiled-{static,shared}`   |
| Fortran local[^19]     | (NYI)          | no               | `ENABLE_FORTRAN=ON`           | `ENABLE_FORTRAN=OFF`                      | `libint_f`                           |

[^15]: Targets for library consumer use. These are available after `find_package(Libint2)` or `add_subdirectory()`.
[^16]: Ensure target found in installation after `find_package(Libint2 COMPONENTS ...)`.
[^17]: Targets in src/lib/libint/CMakeLists.txt.export . Names subject to change. Use namespaced target names in any consuming code.
[^18]: Deprecated legacy aliases. Update any uses to namespaced target.
[^19]: The `libint_f` internal target defines the Fortran interface to Libint2. One must also link to `Libint2::int2` or `Libint2::cxx`. At present, it is not exported, and a namespaced target is not defined.


-----------------------------------------------------------------------------

## Packagers

* Decide if you want the Boost preprocessor headers bundled with Libint (library install includes
  copies of Boost headers) or if they should be a
  build-against-time dependency of the C++11 interface. Withhold (bundle) or supply (dependency)
  Boost detection paths from the library build accordingly. FWIW, Conda bundles.
* Decide if you want the compiled cxx library. something like it is in use in mpqc4


-----------------------------------------------------------------------------

# Platform-Specific Notes

### Linux

### macOS

* Apple `clang++` and [MacPorts](http://www.macports.org/) `g++` (4.8) both work with `-std=c++11` flag
* MacPorts gmp package works fine
* On macOS the default `ar` program lacks support for response files (e.g., https://github.com/evaleev/libint/issues/135 and see https://gitlab.kitware.com/cmake/cmake/issues/16731). Thus you should install the GNU `ar` program (e.g., using HomeBrew: `brew install binutils`) and tell CMake to use it (e.g., add `-DCMAKE_AR=/usr/local/opt/binutils/bin/ar` to the CMake command line).

### Windows

* Several blocking or correctness issues exist; the most thorough list is at .github/workflows/cmake.yml
* A production path is to generate an export tarball with Linux, build static library on Windows, and consume
* Use MPIR package for GMP


-----------------------------------------------------------------------------

# Program-Specific Notes

### mpqc4

* standard libtool configuration:

  ```
  --enable-generic-code --with-max-am=6 --with-opt-am=3 --enable-eri3=0 --enable-eri2=0 --enable-eri3-pure-sh --enable-eri2-pure-sh --enable-fma --disable-1body-property-derivs
  ```

* libtool configuration prior to Jan. 8, 2015:

  ```
  --enable-eri=0 --with-max-am=7 --with-opt-am=4 --disable-unrolling --enable-generic-code --enable-contracted-ints
  ```

### gamess

* standard libtool configuration:

  ```
  --enable-eri=0 --with-max-am=7 --with-opt-am=4 --disable-unrolling --enable-generic-code --enable-contracted-ints --with-cartgauss-ordering=gamess
  ```

### orca

* a libint library (version 2.0.2) is embedded in ORCA
* standard libtool configuration:

  ```
  --enable-eri=2 --enable-eri3=2 --enable-eri2=2 --with-max-am=7 --with-opt-am=4 --with-eri-max-am=7,4,3 --with-eri-opt-am=4,3,2 --disable-unrolling --enable-generic-code --enable-contracted-ints --with-cartgauss-ordering=orca --with-shell-set=orca --enable-eri3-pure-sh --enable-eri2-pure-sh
  ```

### bagel

* standard libtool configuration:

  ```
  --with-max-am=4 --with-eri3-max-am=6 --with-eri2-max-am=6 --enable-eri3=1 -enable-eri=1 --enable-eri2=1 --disable-unrolling --enable-generic-code --enable-contracted-ints --with-cartgauss-ordering=bagel

  ```

* if you want to use spherical Gaussians only add: `--enable-eri3-pure-sh --enable-eri2-pure-sh` (some tests may fail)
* It appears that on a Mac Libint and BAGEL must be either both static or both shared (2/3/2014)

### psi4

* production CMake configuration:

  ```
  -D REQUIRE_CXX_API=ON
  -D LIBINT2_SHGAUSS_ORDERING=gaussian
  -D LIBINT2_CARTGAUSS_ORDERING=standard
  -D LIBINT2_SHELL_SET=standard
  -D ERI3_PURE_SH=OFF
  -D ERI2_PURE_SH=OFF
  -D ENABLE_ERI=2
  -D ENABLE_ERI3=2
  -D ENABLE_ERI2=2
  -D ENABLE_ONEBODY=2
  -D MULTIPOLE_MAX_ORDER=4
  -D WITH_MAX_AM="6;5;4"
  -D WITH_ERI_MAX_AM="5;4;3"
  -D WITH_ERI3_MAX_AM="6;5;4"
  -D WITH_ERI2_MAX_AM="6;5;4"
  ```

* minimal detection:

  ```
  find_package(
    Libint2
    COMPONENTS
      CXX_ho
      gss
      impure_sh
      eri_c4_d0_l3  eri_c3_d0_l4  eri_c2_d0_l4  onebody_d0_l4
      eri_c4_d1_l2  eri_c3_d1_l3  eri_c2_d1_l3  onebody_d1_l3
      eri_c4_d2_l2  eri_c3_d2_l3  eri_c2_d2_l3  onebody_d2_l3
    )
  ```

* see [notes](https://github.com/psi4/psi4/blob/master/external/upstream/libint2/CMakeLists.txt) for details of Psi4/Libint2 configuration

* Psi4 near-future configuration additions:

  ```
    -DDISABLE_ONEBODY_PROPERTY_DERIVS=OFF
    -DENABLE_G12=1
    -DWITH_G12_MAX_AM=4
    -DWITH_G12_OPT_AM=3
  ```


-----------------------------------------------------------------------------

# Miscellaneous Questions

#### Where do I get the source code?

The only way to get the compiler source is from the [Libint source code repository](https://github.com/evaleev/libint) on [GitHub](github.com). You can use a client, like GitHub app or (our favorite) [SourceTree](http://www.sourcetreeapp.com) app from Atlassian. Or from the command line: `git clone https://github.com/evaleev/libint.git`

#### What happened to autoconf?

Version 2.5.0 and older of the exported libint library was buildable using GNU Autoconf and [GNU Make](https://www.gnu.org/software/make/). *As of version 2.6.0 the Autoconf build is deprecated*; the exported libint library should be configured with [CMake](https://cmake.org/) and built with [any CMake-supported generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html), e.g. Ninja and GNU Make.

Version 2.7 and older of the compiler repo was buildable using GNU Autoconf. *As of version 2.8, the Autoconf build is deprecated*; use CMake instead.  TODO 2.8

#### What is the status and importance of SIMD vectorization in Libint?

SIMD vectorization is the crucial contributor to performance of a modern processor core. Libint code can typically hit up to 70% of FLOP peak on a scalar core, hence on a SIMD core divide that number by the vector length (4 for AVX in double precision). The situation is only going to get worse (accelerators already use 8- and 16-wide vector units, and future mainstream processors are likely to use 8-wide units also). Hence if your method spends significant portion of its time computing integrals start rewriting your code now.

Vectorization of Libint is work in progress. However, by switching to AVX we see a factor of 2-2.5 speedup of the integrals kernels compared to scalar performance, thus we are optimistic that it will be possible to attain 50% of peak on AVX hardware. It is clear that significant reorganization of the manner in which integrals are computed and digested is involved, but these costs are unavoidable.

#### What compiler is best?

To obtain peak performance it is *very important* to use the C++ compiler and compiler options that are appropriate for the given platform. It is impossible to provide specific recommendations for specific platforms. The `ENABLE_XHOST` option does allow the compiler to optimize for current architecture. We recommend to use a vendor compiler (e.g., Intel) before trying clang++ and g++. In some situations, however, clang++ and g++ are known to outperform the x86 vendor compiler, so we recommend trying several compilers.


