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

