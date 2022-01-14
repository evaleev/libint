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


-----------------------------------------------------------------------------
## configuring libint compiler

These are the most useful configure options:



###  Which Integrals Classes, Which Derivative Levels

* `ENABLE_ONEBODY` -- Compile with support for up to N-th derivatives of 1-body integrals. Use -1 for OFF. [Default=0]
* `ENABLE_ERI` -- Compile with support for up to N-th derivatives of 4-center electron repulsion integrals. Use -1 for OFF. [Default=0]
* `ENABLE_ERI3` -- Compile with support for up to N-th derivatives of 3-center electron repulsion integrals. Use -1 for OFF. [Default=-1]
* `ENABLE_ERI2` -- Compile with support for up to N-th derivatives of 2-center electron repulsion integrals. Use -1 for OFF. [Default=-1]
* `ENABLE_G12` -- Compile with support for N-th derivatives of MP2-F12 energies with Gaussian factors. Use -1 for OFF. [Default=-1]
* `ENABLE_G12DKH` -- Compile with support for N-th derivatives of DKH-MP2-F12 energies with Gaussian factors. Use -1 for OFF. [Default=-1]


###  How High Angular Momentum

#  example for "semicolon-separated string": `-DENABLE_ERI3=2 -DWITH_ERI3_MAX_AM="5;4;3"`

#  high MAX_AM generating >20k files may require `ulimit -s 65535` for linking library target on Linux to avert "ld: Argument list too long"

* `WITH_MAX_AM` -- Support Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=4]
* `WITH_OPT_AM` -- Optimize maximally for up to angular momentum N (N <= WITH_MAX_AM). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `(WITH_MAX_AM/2)+1`]

* `WITH_ONEBODY_MAX_AM` -- Support 1-body ints for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ONEBODY_OPT_AM` -- Optimize 1-body ints maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI_MAX_AM` -- Support 4-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI_OPT_AM` -- Optimize 4-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI3_MAX_AM` -- Support 3-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI3_OPT_AM` -- Optimize 3-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI2_MAX_AM` -- Support 2-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI2_OPT_AM` -- Optimize 2-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]



### Autotools Update Guide

* `--enable-1body=N` --> `-D ENABLE_ONEBODY`
* `--enable-eri=N` --> `-D ENABLE_ERI=N`
* `--disable-eri` --> `-D ENABLE_ERI=-1`
* `--enable-eri3=N` --> `-D ENABLE_ERI3=N`
* `--enable-eri2=N` --> `-D ENABLE_ERI2=N`

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

* `--enable-shared` --> `-D BUILD_SHARED_LIBS=ON` (standard CMake variable)
* `--enable-static` --> `-D BUILD_SHARED_LIBS=OFF` (standard CMake variable)

* Targets
  * `libint2` --> `Libint2::int2` (internal target name `int-shared`)
  * `libint2_cxx` --> `Libint2::cxx` (internal target name `int-cxx-shared`)


