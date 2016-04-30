#!/bin/sh

set -e

if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
    export OPENMPFLAGS=-fopenmp
else
    # no OpenMP support in clang, will use C++11 threads
    export OPENMPFLAGS=
    export CC=/usr/bin/clang-3.7
    export CXX=/usr/bin/clang++-3.7
fi
export CXXFLAGS="-std=c++11 -Wno-enum-compare $OPENMPFLAGS"
export LDFLAGS=$OPENMPFLAGS
export LIBINT_NUM_THREADS=2

./autogen.sh
./configure CPPFLAGS='-I/usr/include/eigen3' --with-max-am=2,2 --with-eri-max-am=2,1 --with-eri3-max-am=3,2 --enable-eri=1 --enable-eri3=1 --enable-1body=1
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
./configure CPPFLAGS='-I/usr/include/eigen3'
make -j2
make check
