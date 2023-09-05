### Run-Time Compatibility

Functions are provided to check the library configuration and solid harmonics orderings at runtime:

Note: As of v2.8.0 (libtool-based), the configuration_accessor() function will return `(nyi)` by default.
  Packagers are encouraged to patch a generated configuration string into file `configuration.cc` to
  imitate future cmake-based behavior.

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
Configuration: eri_c4_d0_l2;eri_c4_d0_l3;sss;...
```

For the C library, a similar function is available:

```
printf("CMake Configuration (C)  : %s\n", configuration_accessor());
```
```
CMake Configuration (C)  : eri_c4_d0_l2;eri_c4_d0_l3;sss;...
```

If you have a built libint2 library whose history you don't know, a command like this on Linux can provide the same information:

```
strings -n80 /a/random/L2/lying/around/libint2.so
```
```
eri_c2_d0_l2;eri_c2_d0_l3;eri_c2_d1_l2;eri_c3_d0_l2;eri_c3_d0_l3;eri_c3_d1_l2;eri_c4_d0_l2;eri_c4_d1_l2;impure_sh;onebody_d0_l2;onebody_d0_l3;onebody_d1_l2;sss
```

