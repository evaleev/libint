# Docker Travis container notes
This method requires Docker installed on your local machine. This also assumes that you start at the top of the TTG source tree.
1. Create a Travis-CI 'Xenial' docker image: `cd bin; ./docker-travis-build.sh`
2. Run shell in a container using the newly created image: `docker run -it libint-travis-debug bash -l`; if you want to run debugger you want to add `--privileged` flag, as in: `docker run --privileged -it libint-travis-debug bash -l`
3. `cd /home/travis/_build`
4. Configure the job to use the appropriate compiler, compiler version, and debug/release build type:
  * `export BUILD_TYPE=B`, where `B` is `Debug` or `Release`
  * If want to use GNU C++ compiler (gcc):
    * `export GCC_VERSION=VVV` where `VVV` should be the GCC version to be used. The currently valid values are `5`, `6`, `7`, and `8`.
    * `export CXX=g++`
  * If want to use Clang C++ compiler (clang++):
    * `export GCC_VERSION=8`
    * `export CLANG_VERSION=VVV` where `VVV` should be the Clang version to be used. The currently valid value is `8`.
    * `export CXX=clang++`
5. Build and run tests: `./build.sh`
