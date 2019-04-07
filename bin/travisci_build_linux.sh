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
    export CC=/usr/bin/clang-$CLANG_VERSION
    export CXX=/usr/bin/clang++-$CLANG_VERSION
    export FC=/usr/bin/gfortran-$GCC_VERSION
    export EXTRAFLAGS="-stdlib=libc++"
fi
export CXXFLAGS="-std=c++11 -Wno-enum-compare $OPENMPFLAGS $EXTRAFLAGS"
export LDFLAGS=$OPENMPFLAGS
export LIBINT_NUM_THREADS=2

cd ${TRAVIS_BUILD_DIR}
./autogen.sh
cd ${BUILD_PREFIX}
mkdir -p build
cd build
${TRAVIS_BUILD_DIR}/configure CPPFLAGS='-I/usr/include/eigen3' --with-max-am=2,2 --with-eri-max-am=2,2 --with-eri3-max-am=3,2 --enable-eri=1 --enable-eri3=1 --enable-1body=1 --disable-1body-property-derivs --with-multipole-max-order=2
make -j2
make check
cd src/bin/test_eri; ./stdtests.pl; cd ../../..
make export

cd ${BUILD_PREFIX}
mkdir -p export_build
mv build/libint-*.tgz export_build
cd export_build
tar -xvzf libint-*.tgz
rm -f libint-*.tgz
cd libint-*
./configure CPPFLAGS='-I/usr/include/eigen3 -DLIBINT2_DISABLE_BOOST_CONTAINER_SMALL_VECTOR=1' --enable-fortran
make -j2
make check
# build F03 interface
make fortran
fortran/fortran_example
