# Configuring Libint

* Notes
  * Codes "G", "L", or "C" for each option indicate whether it is consumed by the _g_enerator, the _l_ibrary, the library _c_onsumer, or a combination.
    * If your final target is the export tarball, use options that include the letter "G".
    * If you're building a library from an export TARBALL, use options that include the letter "L".
    * For a continuous generator->export->library build, options supplied at the top level will be properly handed off to generator and library build.
  * See [Update Guide](gnu-autotools-update-guide) for new names for old options.


###  Which Integrals Classes, Which Derivative Levels (G)

* `ENABLE_ONEBODY` — G — Compile with support for up to N-th derivatives of 1-body integrals. Use -1 for OFF. [Default=0]
* `ENABLE_ERI` — G — Compile with support for up to N-th derivatives of 4-center electron repulsion integrals. Use -1 for OFF. [Default=0]
* `ENABLE_ERI3` — G — Compile with support for up to N-th derivatives of 3-center electron repulsion integrals. Use -1 for OFF. [Default=-1]
* `ENABLE_ERI2` — G — Compile with support for up to N-th derivatives of 2-center electron repulsion integrals. Use -1 for OFF. [Default=-1]
* `ENABLE_G12` — G — Compile with support for N-th derivatives of MP2-F12 energies with Gaussian factors. Use -1 for OFF. [Default=-1]
* `ENABLE_G12DKH` — G — Compile with support for N-th derivatives of DKH-MP2-F12 energies with Gaussian factors. Use -1 for OFF. [Default=-1]

* `DISABLE_ONEBODY_PROPERTY_DERIVS` — G — Disable geometric derivatives of 1-body property integrals (all but overlap, kinetic, elecpot).
   These derivatives are disabled by default to save compile time. Use OFF to enable.
   Note that the libtool build won't enable this- if forcibly enabled, build_libint balks. [Default=ON]

* `ENABLE_T1G12_SUPPORT` — G — Enable [Ti,G12] integrals when G12 integrals are enabled. Irrelevant when `ENABLE_G12=OFF`. Use OFF to disable. [Default=ON]


###  Which Ordering Conventions (G)

* `LIBINT2_SHGAUSS_ORDERING` — G — Ordering for shells of solid harmonic Gaussians. [Default=standard]
  * `standard` — standard ordering (-l, -l+1 ... l)
  * `gaussian` — the Gaussian ordering (0, 1, -1, 2, -2, ... l, -l)
  See [Solid Harmonic Ordering Scope and History](solid-harmonic-ordering-scope-and-history)
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

####  Solid Harmonic Ordering Scope and History

The use of the `LIBINT2_SHGAUSS_ORDERING` option has changed recently (Dec 2023). See also discussion in the "2023-12-12: 2.8.0" section of [CHANGES](https://github.com/evaleev/libint/blob/master/CHANGES) and in the "cdbb9f3" comment of [engine.h](https://github.com/evaleev/libint/blob/master/include/libint2/engine.h).

Previous to v2.8.0, `LIBINT2_SHGAUSS_ORDERING` was set at generator-build-time for `Operator::sphemultipole` but was re-set-able at library-build-time for other integral classes for both the C and C++ interfaces. Certain macros, `INT_SOLIDHARMINDEX` and `FOR_SOLIDHARM` depended on the library-build-time choice for both the C and C++ interfaces.

Starting at v2.8.0, the generator-build-time construction of `Operator::sphemultipole` has been _fixed at Standard ordering_; thus, the `LIBINT2_SHGAUSS_ORDERING` setting does not influence it. For the C++ interface, solid harmonic ordering for other integrals classes can be toggled at library runtime through `libint2::set_solid_harmonics_ordering`, and the "SOLIDHARM" macros have different forms for accessing either ordering; thus the library build-time `LIBINT2_SHGAUSS_ORDERING` setting is non-constraining. For the C interface, solid harmonic ordering for other integrals classes and the output of "SOLIDHARM" macros _is constrained_ by the `LIBINT2_SHGAUSS_ORDERING` setting. Currently, this setting is fixed at generator-build-time. The build system could be adjusted so that it's re-set-able at library-build-time, but the C interface is discouraged anyways.

Note that options, docs, and CMake components are focused on the C++ interface, and the only remaining constraining influence of the `LIBINT2_SHGAUSS_ORDERING` option -- for the C interface -- may never be acknowledged beyond the previous paragraph.


###  How High Angular Momentum (G)

* Notes
  * example for "semicolon-separated string": `-DENABLE_ERI3=2 -DWITH_ERI3_MAX_AM="5;4;3"`. cmake configuration prints:

    ```
    -- Setting option ENABLE_ERI3: 2
    -- Setting option WITH_ERI3_MAX_AM: 5;4;3
    ```

  * special considerations for high-AM library (L) builds:
    * high MAX_AM generates a large number of source files. If unity builds are disabled, more than
      ~20k files may require `ulimit -s 65535` for linking the library target on Linux to avert
      "ld: Argument list too long".
    * Ninja builds use beyond max threads and can run out of memory, resulting in errorless stops or
      "CMake Error: Generator: execution of make failed". Throttle it to physical threads with
      `export CMAKE_BUILD_PARALLEL_LEVEL=N`.

* `WITH_MAX_AM` — G — Support Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. Specify values greater or equal to `WITH_<class>_MAX_AM`; often mirrors `WITH_ERI3_MAX_AM`. [Default=4]
* `WITH_OPT_AM` — G — Optimize maximally for up to angular momentum N (N <= WITH_MAX_AM). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `(WITH_MAX_AM/2)+1`]

* `MULTIPOLE_MAX_ORDER` — G — Maximum order of spherical multipole integrals. There is no maximum. [Default=4]

* `WITH_ONEBODY_MAX_AM` — G — Support 1-body ints for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ONEBODY_OPT_AM` — G — Optimize 1-body ints maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI_MAX_AM` — G — Support 4-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI_OPT_AM` — G — Optimize 4-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]

* `WITH_ERI3_MAX_AM` — G — Support 3-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. This option controls only the single fitting center; the paired centers use WITH_MAX_AM. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI3_OPT_AM` — G — Optimize 3-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]
* `ERI3_PURE_SH` — G — Assume the 'unpaired' center of 3-center ERIs will be transformed to pure solid harmonics. [Default=OFF]

* `WITH_ERI2_MAX_AM` — G — Support 2-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_ERI2_OPT_AM` — G — Optimize 2-center ERIs maximally for up to angular momentum N (N <= max-am). Can specify values for each derivative level as a semicolon-separated string. [Default=-1 -> `WITH_OPT_AM`]
* `ERI2_PURE_SH` — G — Assume the 2-center ERIs will be transformed to pure solid harmonics. [Default=OFF]

* `WITH_G12_MAX_AM` — G — Support integrals for G12 methods of angular momentum up to N. No specification with per-derivative list. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_G12_OPT_AM` — G — Optimize G12 integrals for up to angular momentum N (N <= max-am). No specification with per-derivative list. [Default=-1 `WITH_OPT_AM`]

* `WITH_G12DKH_MAX_AM` — G — Support integrals for relativistic G12 methods of angular momentum up to N. No specification with per-derivative list. [Default=-1 -> `WITH_MAX_AM`]
* `WITH_G12DKH_OPT_AM` — G — Optimize G12DKH integrals for up to angular momentum N (N <= max-am). No specification with per-derivative list. [Default=-1 `WITH_OPT_AM`]


### Compilers and Flags (G L) (TARBALL)

### Miscellaneous (G L)

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

* `LIBINT_USER_DEFINED_REAL_INCLUDES` — L — Additional #includes necessary to use the real type. [Defaults=none]
* `LIBINT_CONTRACTED_INTS` — G — Turn on support for contracted integrals. [Default=ON]
* `LIBINT_ERI_STRATEGY` — G — Compute ERIs using the following strategy (experts only). (0 for OS, 1 for HGP, 2 for HL). [Default=1]
* `LIBINT_USE_COMPOSITE_EVALUATORS` — G — Libint will use composite evaluators (i.e. every evaluator will compute one integral type only). [Default=ON]
* `LIBINT_SINGLE_EVALTYPE` — G — Generate single evaluator type (i.e. all tasks use the same evaluator). OFF is NYI [Default=ON]
* `LIBINT_ENABLE_UNROLLING` — G — Unroll shell sets into integrals (will unroll shell sets larger than N) (0 for never, N for N, 1000000000 for always). [Default=100]
* `LIBINT_ALIGN_SIZE` — G — If posix_memalign is available, this will specify alignment of Libint data, in units of sizeof(LIBINT2_REALTYPE). Default is to use built-in heuristics: system-determined for vectorization off (default) or veclen * sizeof(LIBINT2_REALTYPE) for vectorization on. (experts only). [Default=0]
* `LIBINT_GENERATE_FMA` — G — Generate FMA (fused multiply-add) instructions (to benefit must have FMA-capable hardware and compiler). [Default=OFF]
* `LIBINT_ENABLE_GENERIC_CODE` — G — Use manually-written generic code. [Default=OFF]
* `LIBINT_API_PREFIX` — G — Prepend this string to every name in the library API (except for the types). [Default=OFF]
* `LIBINT_VECTOR_LENGTH` — G — Compute integrals in vectors of length N. [Default=OFF]
* `LIBINT_VECTOR_METHOD` — G — Specifies how to vectorize integrals. Irrelevant when `LIBINT_VECTOR_LENGTH=OFF. Allowed values are 'block' and 'line'.  [Default=block]
* `LIBINT_ACCUM_INTS` — G — Accumulate integrals to the buffer, rather than copy (OFF for copy, ON for accum). [Default=OFF]
* `LIBINT_FLOP_COUNT` — G — Support (approximate) FLOP counting by the library. (Generated code will require C++11!). [Default=OFF]
* `LIBINT_PROFILE` — G — Turn on profiling instrumentation of the library. (Generated code will require C++11!). [Default=OFF]


# GNU Autotools Update Guide

* Notes
  * Multiple option names can be from any long-lived branch but usually libtool+cmake --> final cmake+cmake.

* `--enable-1body=N` --> `-D ENABLE_ONEBODY=N`
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


* Defunct as non-CMake-like: `--build`, `--host`, `--target`,

* `--with-api-prefix=pfx` --> `-D LIBINT_API_PREFIX=pfx`
* `--enable-unrolling=yes` --> `-D LIBINT_ENABLE_UNROLLING=1000000000`
* `--enable-unrolling=no` --> `-D LIBINT_ENABLE_UNROLLING=0`
* `--enable-unrolling=S` --> `-D LIBINT_ENABLE_UNROLLING=S`
* `--enable-generic-code` --> `-D LIBINT_ENABLE_GENERIC_CODE=ON`
* `--with-vector-length=N` --> `-D LIBINT_VECTOR_LENGTH=N`
* `--with-vector-method=choice` --> `-D LIBINT_VECTOR_METHOD=choice`
* `--with-align-size=N` --> `-D LIBINT_ALIGN_SIZE=N` (G) (`-D LIBINT2_ALIGN_SIZE=N` for L)
* `--enable-fma` --> `-D LIBINT_GENERATE_FMA=ON`
* `--enable-accum-ints` --> `-D LIBINT_ACCUM_INTS=ON`
* `--enable-flop-counter` -> `-D LIBINT_FLOP_COUNT=ON`
* `--enable-profile` --> `-D LIBINT_PROFILE=ON`
* `--disable-contracted-ints` --> `-D LIBINT_CONTRACTED_INTS=OFF`
* `--disable-single-evaltype` --> `-D LIBINT_SINGLE_EVALTYPE=OFF` (NYI)
* `--enable-composite-evaluators` --> `-D LIBINT_USE_COMPOSITE_EVALUATORS=ON`
* `--disable-composite-evaluators` --> `-D LIBINT_USE_COMPOSITE_EVALUATORS=OFF`
* `--with-eri-strategy=OS` --> `-D LIBINT_ERI_STRATEGY=0`
* `--with-eri-strategy=HL` --> `-D LIBINT_ERI_STRATEGY=2`
* `--with-real-type=type` --> `-D LIBINT2_REALTYPE=type`
* `--with-real-type-inclues=inc` --> `-D LIBINT_USER_DEFINED_REAL_INCLUDES="#include <stdio.h>"`





### Run-Time Compatibility

Functions are provided to check the library configuration and solid harmonics orderings at runtime:

Note: As of v2.8.0 (libtool-based), the configuration_accessor() function will return `(nyi)` by default.
  Packagers are encouraged to patch a generated configuration string into file `configuration.cc.cmake.in` to
  imitate future cmake-based behavior. See sample patch below. The string can be generated by editing
  and running `export/cmake/configuration-gen.py`.
  Also patch MAX_AM_ERI in `CMakeLists.txt` (of export tarball; `export/cmake/CMakeLists.txt.export` in repo src).

```
libint2::initialize();
printf("SHGShell: %d\n", libint2::solid_harmonics_ordering());
libint2::set_solid_harmonics_ordering(libint2::SHGShellOrdering_Gaussian);
// note that toggling solid_harmonics_ordering works fine for printing here, but
//    it's recc. to set *before* Engines are created, so before calling initialize()
printf("SHGShell: %d\n", libint2::solid_harmonics_ordering());
// if patched as described above
printf("Configuration: %s\n", libint2::configuration_accessor().c_str());
printf("Supports: dddd=%d mmmm=%d\n", libint2::supports("eri_dddd_d0"), libint2::supports("eri_mmmm_d0"));
libint2::finalize();
```
```
SHGShell: 1
SHGShell: 2
Configuration: eri_dddd_d0_l2;eri_ffff_d0;ss;...
Supports: dddd=1 mmmm=0
```

For the C library, a similar function is available:

```
printf("CMake Configuration (C)  : %s\n", configuration_accessor());
```
```
CMake Configuration (C)  : eri_dddd_d0;eri_ffff_d0;ss;...
```

If you have a built libint2 library whose history you don't know, a command like this on Linux can provide the same information:

```
strings -n80 /a/random/L2/lying/around/libint2.so
```
```
ss;onebody_ii_d0;onebody_hh_d0;eri_iiI_d0;eri_iii_d0;eri_II_d0;eri_ii_d0;eri_hhhh_d0;eri_hhH_d0;eri_hhh_d0;eri_HH_d0;eri_hh_d0;eri_gggg_d1;eri_dddd_d1
```

A patch like the following is suitable for an export tarball generated from the next following.
[See guide](#configuration-codes) for decoding the configuration components.

```
--- src/configuration.cc.cmake.in   2023-09-05 09:13:50.000000000 -0400
+++ src/configuration.cc.cmake.in_basic 2023-09-05 23:41:00.444396591 -0400
@@ -24,6 +24,6 @@
    @return the semicolon-separated strings from CMake components */
 const char * configuration_accessor() {
     //return "@Libint2_CONFIG_COMPONENTS@";
-    return "(nyi)";
+    return "ss;multipole_nn_d0;multipole_mm_d0;multipole_ll_d0;multipole_kk_d0;multipole_ii_d0;multipole_hh_d0;multipole_gg_d0;multipole_ff_d0;multipole_dd_d0;onebody_ii_d0;onebody_hh_d0;onebody_gg_d0;onebody_ff_d0;onebody_dd_d0;onebody_hh_d1;onebody_gg_d1;onebody_ff_d1;onebody_dd_d1;onebody_gg_d2;onebody_ff_d2;onebody_dd_d2;eri_hhhh_d0;eri_gggg_d0;eri_ffff_d0;eri_dddd_d0;eri_gggg_d1;eri_ffff_d1;eri_dddd_d1;eri_iiI_d0;eri_hhI_d0;eri_hhH_d0;eri_ggI_d0;eri_ggH_d0;eri_ggG_d0;eri_ffI_d0;eri_ffH_d0;eri_ffG_d0;eri_ffF_d0;eri_ddI_d0;eri_ddH_d0;eri_ddG_d0;eri_ddF_d0;eri_ddD_d0;eri_hhH_d1;eri_ggH_d1;eri_ggG_d1;eri_ffH_d1;eri_ffG_d1;eri_ffF_d1;eri_ddH_d1;eri_ddG_d1;eri_ddF_d1;eri_ddD_d1;eri_iii_d0;eri_hhi_d0;eri_hhh_d0;eri_ggi_d0;eri_ggh_d0;eri_ggg_d0;eri_ffi_d0;eri_ffh_d0;eri_ffg_d0;eri_fff_d0;eri_ddi_d0;eri_ddh_d0;eri_ddg_d0;eri_ddf_d0;eri_ddd_d0;eri_hhh_d1;eri_ggh_d1;eri_ggg_d1;eri_ffh_d1;eri_ffg_d1;eri_fff_d1;eri_ddh_d1;eri_ddg_d1;eri_ddf_d1;eri_ddd_d1;eri_II_d0;eri_HH_d0;eri_GG_d0;eri_FF_d0;eri_DD_d0;eri_HH_d1;eri_GG_d1;eri_FF_d1;eri_DD_d1;eri_ii_d0;eri_hh_d0;eri_gg_d0;eri_ff_d0;eri_dd_d0;eri_hh_d1;eri_gg_d1;eri_ff_d1;eri_dd_d1;g12_gggg_d0;g12_ffff_d0;g12_dddd_d0;g12_gggg_d1;g12_ffff_d1;g12_dddd_d1";
}
```
```
./configure \
    --enable-eri=1 \
    --enable-eri3=1 \
    --enable-eri2=1 \
    --enable-1body=2 \
    --enable-g12=1 \
    --disable-1body-property-derivs \
    --with-multipole-max-order=10 \
    --with-g12-max-am=4 \
    --with-eri-max-am=5,4 \
    --with-eri3-max-am=6,5 \
    --with-eri2-max-am=6,5 \
    --with-max-am=6,5
```


#### Configuration Codes

Eventually, these will be CMake Components, too.

```
   multipole_hh_dD - library includes spherical multipole integrals with max angular momentum up to
                     "h" (h=spdfghikl...; s,p not enumerated) and derivative order "D" (D=0,1,2,...).
                     For example, the presence of "multipole_ii_d0" means mpole ints are available for L=6.
   onebody_hh_dD   - library includes 1-body integrals with max angular momentum up to "h"
                     (h=spdfghikl...; s,p not enumerated) and derivative order "D" (D=0,1,2,...).
                     For example, the presence of "onebody_ii_d1" means onebody gradient ints are
                     available for L=6.
   eri_hhhh_dD     - library includes 2-body integrals with 4 centers and max angular momentum up to
                     "h" (h=spdfghikl...; s,p not enumerated) and derivative order "D" (D=0,1,2,...).
                     For example, the presence of "eri_ffff_d1" means 4-center gradient ints are
                     available for L=3. That is, the library was configured with at least
                     '-D ENABLE_ERI=1 -D WITH_ERI_MAX_AM="?;>=3"'.
   eri_hhL_dD      - library includes 2-body integrals with 3 centers and max angular momentum up to
   eri_hhl_dD        Cartesian "h" for the two paired centers and Cartesian "l" or solid harmonics "L"
                     for the unpaired/fitting center, (h/l=spdfghikl..., L=SPDFGHIKL...; l>=h
                     enumerated; s,p,S,P not enumerated) and derivative order "D" (D=0,1,2,...). The
                     "eri_hhL_dD" component is always available when 3-center ints are present. When pure
                     solid harmonics are assumed for 3-center ints, "eri_hhl_dD" will *not be available*.
                     For example, the presence of "eri_ffG_d0" means 3-center energy ints are
                     available for L=3 (paired centers) and L=4 (fitting center). That is, the library
                     was configured with at least "-D ENABLE_ERI3=0 -D WITH_MAX_AM=3 -D WITH_ERI3_MAX_AM=4".
                     The presence of "eri_ffg_d0" means the library configuration did not additionally
                     include "-D ERI3_PURE_SH=ON".
   eri_HH_dD       - library includes 2-body integrals with 2 centers and max angular momentum up to
   eri_hh_dD         Cartesian "h" or solid harmonics "H", (h=spdfghikl..., H=SPDFGHIKL...; s,p,S,P not
                     enumerated) and derivative order "D" (D=0,1,2,...). The "eri_HH_dD" component is
                     always available when 2-center ints are present. When pure solid harmonics are
                     assumed for 2-center ints, "eri_hh_dD" will *not be available*.
                     For example, the presence of "eri_FF_d2" means 2-center Hessian ints are
                     available for L=3. That is, the library was configured with at least
                     '-D ENABLE_ERI2=2 -D WITH_ERI2_MAX_AM="?;?;>=3"'. The presence of "eri_ff_d2" means the
                     library configuration did not additionally include "-D ERI2_PURE_SH=ON".
   g12_hhhh_dD     - library includes F12 integrals with Gaussian factors and max angular momentum up to
                     "h" (h=spdfghikl...; s,p not enumerated) and derivative order "D" (D=0,1,2,...).
                     For example, the presence of "g12_iiii_d2" means g12 Hessian ints are available for L=6.

                                        cart       shell_set   used_by
                                        --------   ---------   -------
   ss - library integrals use ordering  standard + standard  = mpqc4, cp2k, psi4 (psi4 requires runtime-setting of solid harmonic ordering to Gaussian)
   so - library integrals use ordering           + orca
   is - library integrals use ordering  intv3    + standard  = mpqc3
   io - library integrals use ordering           + orca
   gs - library integrals use ordering  gamess   + standard  = gamess
   go - library integrals use ordering           + orca
   os - library integrals use ordering  orca     + standard
   oo - library integrals use ordering           + orca      = orca
   bs - library integrals use ordering  bagel    + standard  = bagel
   bo - library integrals use ordering           + orca
```

### Interfacing

Eventually (approximately 2.9.0 CMake-based), additional functions will be available to retrive Libint version, commit, and literature citation. Below are outputs at the libtool stage.

```
auto Mmp = libint2::libint_version();
printf("Version: Numeric=%s Sortable=%s Commit=%s\n", libint2::libint_version_string(false).c_str(), libint2::libint_version_string(true).c_str(), libint2::libint_commit().c_str());
printf("Version: Major=%d minor=%d patch=%d\n", std::get<0>(Mmp), std::get<1>(Mmp), std::get<2>(Mmp));
printf("Citation: DOI=%s Ref=%s\n", libint2::libint_reference_doi().c_str(), libint2::libint_reference().c_str());
printf("Citation: BibTex=%s\n", libint2::libint_bibtex().c_str());
```
```
Version: Numeric=2.8.0 Sortable= Commit=
Version: Major=2 minor=8 patch=0
Citation: DOI= Ref=Libint: , Version  Edward F. Valeev, http://libint.valeyev.net/
Citation: BibTex=@Misc{Libint2,
  author = {E.~F.~Valeev},
  title = {\textsc{Libint}: },
  howpublished = {http://libint.valeyev.net/},
  note = {version },
  year = {}
}
```
