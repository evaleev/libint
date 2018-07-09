#!/bin/sh

set -e

if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
    export FC=/usr/bin/gfortran-$GCC_VERSION
    export OPENMPFLAGS=-fopenmp
    export EXTRAFLAGS=
else
    # no OpenMP support in clang, will use C++11 threads
    export OPENMPFLAGS=
    export CC=/usr/bin/clang-5.0
    export CXX=/usr/bin/clang++-5.0
    export FC=/usr/bin/gfortran-$GCC_VERSION
    # Boost 1.55 is too old, override Boost.PP detection of variadic macro support
    export EXTRAFLAGS="-DBOOST_PP_VARIADICS=1"
fi
export CXXFLAGS="-std=c++11 -Wno-enum-compare $OPENMPFLAGS $EXTRAFLAGS"
export LDFLAGS=$OPENMPFLAGS
export LIBINT_NUM_THREADS=2

./autogen.sh
./configure CPPFLAGS='-I/usr/include/eigen3' --with-max-am=2,2 --with-eri-max-am=2,2 --with-eri3-max-am=3,2 --enable-eri=1 --enable-eri3=1 --enable-1body=1 --disable-1body-property-derivs --with-multipole-max-order=2
make -j2
make check
cd src/bin/test_eri; ./stdtests.pl; cd ../../..

make export
mkdir export_build
mv libint-*.tgz export_build
cd export_build
tar -xvzf libint-*.tgz
rm -f libint-*.tgz
cd libint-*
./configure CPPFLAGS='-I/usr/include/eigen3' --enable-fortran
make -j2
make check
# build F03 interface
make fortran
fortran/fortran_example
