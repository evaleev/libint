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

* `WITH_MAX_AM` — G — Support Gaussians of angular momentum up to N. If ERI3 ints are enabled, specifing values for each derivative level as a semicolon-separated string also controls the AM of the paired centers. [Default=4]
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


# GNU Autotools Update Guide

* Notes
  * Multiple option names can be from any long-lived branch but usually libtool+cmake --> final cmake+cmake.

* `--enable-1body=N` --> `-D ENABLE_ONEBODY=N`
* `--enable-eri=N` --> `-D ENABLE_ERI=N`
* `--disable-eri` --> `-D ENABLE_ERI=-1`
* `--enable-eri3=N` --> `-D ENABLE_ERI3=N`
* `--enable-eri2=N` --> `-D ENABLE_ERI2=N`

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
                     "h" (h=spdfghikl...) and derivative order "D" (D=0,1,2,...).
                     For example, the presence of "multipole_ii_d0" means mpole ints are available for L=6.
   onebody_hh_dD   - library includes 1-body integrals with max angular momentum up to "h"
                     (h=spdfghikl...) and derivative order "D" (D=0,1,2,...).
                     For example, the presence of "onebody_ii_d1" means onebody gradient ints are
                     available for L=6.
   eri_hhhh_dD     - library includes 2-body integrals with 4 centers and max angular momentum up to
                     "h" (h=spdfghikl...) and derivative order "D" (D=0,1,2,...).
                     For example, the presence of "eri_ffff_d1" means 4-center gradient ints are
                     available for L=3. That is, the library was configured with at least
                     '-D ENABLE_ERI=1 -D WITH_ERI_MAX_AM="?;>=3"'.
   eri_hhL_dD      - library includes 2-body integrals with 3 centers and max angular momentum up to
   eri_hhl_dD        Cartesian "h" for the two paired centers and Cartesian "l" or solid harmonics "L"
                     for the unpaired/fitting center, (h/l=spdfghikl..., L=SPDFGHIKL...; l>=h
                     enumerated) and derivative order "D" (D=0,1,2,...). The
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
                     "h" (h=spdfghikl...) and derivative order "D" (D=0,1,2,...).
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

Eventually (approximately 2.10.0 CMake-based), additional functions will be available to retrive Libint version, commit, and literature citation. Below are outputs at the libtool stage.

```
auto Mmp = libint2::libint_version();
printf("Version: Numeric=%s Sortable=%s Commit=%s\n", libint2::libint_version_string(false).c_str(), libint2::libint_version_string(true).c_str(), libint2::libint_commit().c_str());
printf("Version: Major=%d minor=%d patch=%d\n", std::get<0>(Mmp), std::get<1>(Mmp), std::get<2>(Mmp));
printf("Citation: DOI=%s Ref=%s\n", libint2::libint_reference_doi().c_str(), libint2::libint_reference().c_str());
printf("Citation: BibTeX=%s\n", libint2::libint_bibtex().c_str());
```
```
Version: Numeric=2.8.0 Sortable= Commit=
Version: Major=2 minor=8 patch=0
Citation: DOI= Ref=Libint: , Version  Edward F. Valeev, http://libint.valeyev.net/
Citation: BibTeX=@Misc{Libint2,
  author = {E.~F.~Valeev},
  title = {\textsc{Libint}: },
  howpublished = {http://libint.valeyev.net/},
  note = {version },
  year = {}
}
```
