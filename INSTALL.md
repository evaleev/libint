### Run-Time Compatibility

Functions are provided to check the library configuration and solid harmonics orderings at runtime:

Note: As of v2.8.0 (libtool-based), the configuration_accessor() function will return `(nyi)` by default.
  Packagers are encouraged to patch a generated configuration string into file `configuration.cc` to
  imitate future cmake-based behavior. See sample patch below.

```
libint2::initialize();
printf("SHGShell: %d\n", libint2::solid_harmonics_ordering());
libint2::set_solid_harmonics_ordering(libint2::SHGShellOrdering_Gaussian);
printf("SHGShell: %d\n", libint2::solid_harmonics_ordering());
printf("Configuration: %s\n", libint2::configuration_accessor().c_str());
libint2::finalize();
```
```
SHGShell: 1
SHGShell: 2
Configuration: eri_c4_d0_l2;eri_c4_d0_l3;ss;...
```

For the C library, a similar function is available:

```
printf("CMake Configuration (C)  : %s\n", configuration_accessor());
```
```
CMake Configuration (C)  : eri_c4_d0_l2;eri_c4_d0_l3;ss;...
```

If you have a built libint2 library whose history you don't know, a command like this on Linux can provide the same information:

```
strings -n80 /a/random/L2/lying/around/libint2.so
```
```
eri_c2_d0_l2;eri_c2_d0_l3;eri_c2_d1_l2;eri_c3_d0_l2;eri_c3_d0_l3;eri_c3_d1_l2;eri_c4_d0_l2;eri_c4_d1_l2;impure_sh;onebody_d0_l2;onebody_d0_l3;onebody_d1_l2;ss
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
+    return "eri_c2_d0_l2;eri_c2_d0_l3;eri_c2_d0_l4;eri_c2_d0_l5;eri_c2_d0_l6;eri_c2_d1_l2;eri_c2_d1_l3;eri_c2_d1_l4;eri_c2_d1_l5;eri_c3_d0_l2;eri_c3_d0_l3;eri_c3_d0_l4;eri_c3_d0_l5;eri_c3_d0_l6;eri_c3_d1_l2;eri_c3_d1_l3;eri_c3_d1_l4;eri_c3_d1_l5;eri_c4_d0_l2;eri_c4_d0_l3;eri_c4_d0_l4;eri_c4_d0_l5;eri_c4_d1_l2;eri_c4_d1_l3;eri_c4_d1_l4;g12_d0_l2;g12_d0_l3;g12_d0_l4;g12_d1_l2;g12_d1_l3;g12_d1_l4;impure_sh;onebody_d0_l2;onebody_d0_l3;onebody_d0_l4;onebody_d0_l5;onebody_d0_l6;onebody_d1_l2;onebody_d1_l3;onebody_d1_l4;onebody_d1_l5;onebody_d2_l2;onebody_d2_l3;onebody_d2_l4;ss";
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
    --with-multipole-max-order=4 \
    --with-g12-max-am=4 \
    --with-eri-max-am=5,4 \
    --with-eri3-max-am=6,5 \
    --with-eri2-max-am=6,5 \
    --with-max-am=6,5
```

#### Configuration Codes

Evenually, these will be CMake Components, too.

```
   onebody_dD_lL - library includes 1-body integrals with derivative order D (D=0,1,2,...) and max angular momentum up to L (L=2,3,4,...)
   eri_cC_dD_lL  - library includes 2-body integrals with C (C=2,3,4) centers, derivative order D (D=0,1,2,...), and max angular momentum up to L (L=2,3,4,...)
   g12_dD-lL - library includes F12 integrals with Gaussian factors with derivative order D and max angular momentum up to L

   impure_sh - library doesn't assume that 2- and 3-center integrals involve pure solid harmonics

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
