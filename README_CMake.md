# Libint

### History

This is the Libint project (http://evaleev.github.io/libint/) by
Prof. Edward F. Valeev (@evaleev) with early roots by Prof. Justin T.
Fermann.

Libint1 has source available on sourceforge 
(https://sourceforge.net/projects/libint/files/v1-releases/) or in this repository as a branch on GitHub
(https://github.com/evaleev/libint/tree/v1) and primarily
builds with `make`. Libint1 separates the build process for derivative integrals.

Libint2 has source available on GitHub (https://github.com/evaleev/libint) and,
as distributed, builds with `make`. Libint2 integrates the build process for
derivative integrals.

Libint1 and Libderiv1 have been in the *ab initio* quantum chemistry package Psi4
(http://psicode.org/, https://github.com/psi4/psi4) since 2009. Internal to Psi4, it
has, since about 2014, built with `cmake`, as designed by @andysim.
@ryanmrichard and @loriab have reworked the CMake build and extracted the
project until suitable for `ExternalProject_Add`. As of version 1.2.0,
Psi4's CMake build system for libint1 and libderiv1 source has been
ported back to this main Libint repository as an alternate build system.

#### Building

```bash
cmake -H. -Bobjdir \
 -DMAX_AM_ERI=4
cd objdir && make
make install 
```

The primary CMake option is `MAX_AM_ERI` to control the maximum angular
momentum for integrals. This is a Psi4 quantity slightly different from
those used internally by Libint and found in the installed header files:

`MAX_AM_ERI` | `LIBINT_MAX_AM` | `LIBINT_OPT_AM` | `LIBDERIV_MAX_AM1` | `LIBDERIV_MAX_AM12`
------------ | --------------- | --------------- | ------------------ | -------------------
3 | 4 | 4 | 3 | 2
4 | 5 | 3 | 4 | 3
**5** | **6** | **3** | **5** | **4**
6 | 7 | 4 | 6 | 5
7 | 8 | 4 |
8 | 9 | 5 |

For orientation on an atom such as Neon, the default **5** gets you conventional cc-pV5Z for energies, cc-pVQZ for gradients, cc-pVTZ for frequencies and density-fitted cc-pVQZ for energies, cc-pVTZ for gradients, cc-pVDZ for frequencies.

The build is also responsive to 

* static/shared toggle `BUILD_SHARED_LIBS`
* the install location `CMAKE_INSTALL_PREFIX`
* of course, `CMAKE_C_COMPILER`, `CMAKE_CXX_COMPILER`, `CMAKE_C_FLAGS`, and `CMAKE_CXX_FLAGS`

See [CMakeLists.txt](CMakeLists.txt) for options details and additional options. All these build options should be passed as `cmake -DOPTION`.

#### Detecting

This project installs with `LibintConfig.cmake` and `LibintConfigVersion.cmake` files suitable for use with CMake [`find_package()`](https://cmake.org/cmake/help/v3.2/command/find_package.html) in `CONFIG` mode.

* `find_package(Libint)` - find any Libint libraries and headers
* `find_package(Libint 1.1.6 EXACT CONFIG REQUIRED COMPONENTS shared 6)` - find Libint exactly version 1.1.6 built with shared libraries and `MAX_AM_ERI` >= 6 or die trying

See [LibintConfig.cmake.in](cmake/LibintConfig.cmake.in) for details of how to detect the Config file and what CMake variables and targets are exported to your project.

#### Using

After `find_package(Libint ...)`,

* test if package found with `if(${Libint_FOUND})` or `if(TARGET Libint::libint)`
* link to library (establishes dependency), including header and definitions configuration with `target_link_libraries(mytarget Libint::libint)`
* include header files using `target_include_directories(mytarget PRIVATE $<TARGET_PROPERTY:Libint::libint,INTERFACE_INCLUDE_DIRECTORIES>)`
* compile target applying `-DUSING_Libint;-DMAX_AM_ERI=N` definition using `target_compile_definitions(mytarget PRIVATE $<TARGET_PROPERTY:Libint::libint,INTERFACE_COMPILE_DEFINITIONS>)`
