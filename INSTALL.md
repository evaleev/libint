### Run-Time Compatibility

Functions are provided to check the library configuration and solid harmonics orderings at runtime:

Note: As of v2.8.0 (libtool-based), the configuration_accessor() function will return `(nyi)` by default.
  Packagers are encouraged to patch a generated configuration string into file `configuration.cc.cmake.in` to
  imitate future cmake-based behavior. See sample patch below.

```
libint2::initialize();
printf("SHGShell: %d\n", libint2::solid_harmonics_ordering());
libint2::set_solid_harmonics_ordering(libint2::SHGShellOrdering_Gaussian);
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
ss;onebody_i_d0;onebody_h_d0;eri_iiI_d0;eri_iii_d0;eri_II_d0;eri_ii_d0;eri_hhhh_d0;eri_hhH_d0;eri_hhh_d0;eri_HH_d0;eri_hh_d0;eri_gggg_d1;eri_dddd_d1
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
+    return "ss;multipole_n;multipole_m;multipole_l;multipole_k;multipole_i;multipole_h;multipole_g;multipole_f;multipole_d;onebody_i_d0;onebody_h_d0;onebody_g_d0;onebody_f_d0;onebody_d_d0;onebody_h_d1;onebody_g_d1;onebody_f_d1;onebody_d_d1;onebody_g_d2;onebody_f_d2;onebody_d_d2;eri_hhhh_d0;eri_gggg_d0;eri_ffff_d0;eri_dddd_d0;eri_gggg_d1;eri_ffff_d1;eri_dddd_d1;eri_iiI_d0;eri_hhI_d0;eri_hhH_d0;eri_ggI_d0;eri_ggH_d0;eri_ggG_d0;eri_ffI_d0;eri_ffH_d0;eri_ffG_d0;eri_ffF_d0;eri_ddI_d0;eri_ddH_d0;eri_ddG_d0;eri_ddF_d0;eri_ddD_d0;eri_hhH_d1;eri_ggH_d1;eri_ggG_d1;eri_ffH_d1;eri_ffG_d1;eri_ffF_d1;eri_ddH_d1;eri_ddG_d1;eri_ddF_d1;eri_ddD_d1;eri_iii_d0;eri_hhi_d0;eri_hhh_d0;eri_ggi_d0;eri_ggh_d0;eri_ggg_d0;eri_ffi_d0;eri_ffh_d0;eri_ffg_d0;eri_fff_d0;eri_ddi_d0;eri_ddh_d0;eri_ddg_d0;eri_ddf_d0;eri_ddd_d0;eri_hhh_d1;eri_ggh_d1;eri_ggg_d1;eri_ffh_d1;eri_ffg_d1;eri_fff_d1;eri_ddh_d1;eri_ddg_d1;eri_ddf_d1;eri_ddd_d1;eri_II_d0;eri_HH_d0;eri_GG_d0;eri_FF_d0;eri_DD_d0;eri_HH_d1;eri_GG_d1;eri_FF_d1;eri_DD_d1;eri_ii_d0;eri_hh_d0;eri_gg_d0;eri_ff_d0;eri_dd_d0;eri_hh_d1;eri_gg_d1;eri_ff_d1;eri_dd_d1;g12_g_d0;g12_f_d0;g12_d_d0;g12_g_d1;g12_f_d1;g12_d_d1";
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

Evenually, these will be CMake Components, too.

```
   multipole_h  - library includes spherical multipole integrals with max angular momentum up to
                  "h" (h=(sp)dfghikl...; s,p not enumerated).
                  For example, the presence of "multipole_i" means mpole ints are available for L=6.
   onebody_h_dD - library includes 1-body integrals with max angular momentum up to "h"
                  (h=(sp)dfghikl...; s,p not enumerated) and derivative order "D" (D=0,1,2,...).
                  For example, the presence of "onebody_i_d1" means onebody gradient ints are
                  available for L=6.
   eri_hhhh_dD  - library includes 2-body integrals with 4 centers and max angular momentum up to
                  "h" (h=(sp)dfghikl...; s,p not enumerated) and derivative order "D" (D=0,1,2,...).
                  For example, the presence of "eri_ffff_d1" means 4-center gradient ints are
                  available for L=3. That is, the library was configured with at least
                  "--enable-eri=1 --with-eri-max-am=?,>=3".
   eri_hhL_dD   - library includes 2-body integrals with 3 centers and max angular momentum up to
   eri_hhl_dD     Cartesian "h" for the two paired centers and Cartesian "l" or solid harmonics "L"
                  for the unpaired/fitting center, (h/l=(sp)dfghikl..., L=(SP)DFGHIKL...; l>=h
                  enumerated; s,p not enumerated) and derivative order "D" (D=0,1,2,...). The
                  "eri_hhL_dD" component is always available when 3-center ints are present. When pure
                  solid harmonics are assumed for 3-center ints, "eri_hhl_dD" will *not be available*.
                  For example, the presence of "eri_ffG_d0" means 3-center energy ints are
                  available for L=3 (paired centers) and L=4 (fitting center). That is, the library
                  was configured with at least "--enable-eri3=0 --with-max-am=3 --with-eri3-max-am=4".
                  The presence of "eri_ffg_d0" means the library configuration did not additionally
                  include "--enable-eri3-pure-sh[=yes]".
   eri_HH_dD    - library includes 2-body integrals with 2 centers and max angular momentum up to
   eri_hh_dD      Cartesian "h" or solid harmonics "H", (h=(sp)dfghikl..., H=(SP)DFGHIKL...; s,p not
                  enumerated) and derivative order "D" (D=0,1,2,...). The "eri_HH_dD" component is
                  always available when 2-center ints are present. When pure solid harmonics are
                  assumed for 2-center ints, "eri_hh_dD" will *not be available*.
                  For example, the presence of "eri_FF_d2" means 2-center Hessian ints are
                  available for L=3. That is, the library was configured with at least
                  "--enable-eri2=2 --with-eri2-max-am=?,?,>=3". The presence of "eri_ff_d2" means the
                  library configuration did not additionally include "--enable-eri2-pure-sh[=yes]".
   g12_h_dD     - library includes F12 integrals with Gaussian factors max angular momentum up to
                  "h" (h=(sp)dfghikl...; s,p not enumerated) and derivative order "D" (D=0,1,2,...).
                  For example, the presence of "g12_i_d2" means g12 Hessian ints are available for L=6.

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
