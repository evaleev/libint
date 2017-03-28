#!/bin/sh

set -e

VERSION=1.63.0
VERSION_="$(echo ${VERSION} | sed 's/\./_/g')"
TARBALL=boost_${VERSION_}.tar.bz2

#### download
if ! test -f ${TARBALL}; then
  BOOST_URL=https://sourceforge.net/projects/boost/files/boost/${VERSION}/${TARBALL}/download
  curl -o ${TARBALL} -L ${BOOST_URL} 
fi

#### unpack
if ! test -d boost_${VERSION_}; then
  tar -xvjf ${TARBALL}
fi

#### build bcp
cd boost_${VERSION_}
./bootstrap.sh && ./b2 tools/bcp

#### extract boost/preprocessor.hpp and dependencies
mkdir tmp
bin.v2/tools/bcp/darwin-4.2.1/release/link-static/bcp --unix-lines boost/preprocessor.hpp tmp

#### specialize to libint2 by prepending all includes with libint2
find tmp/boost -type f -print0 | xargs -0 sed -i "" 's/ <boost\// <libint2\/boost\//g'
cd tmp
tar -cvzf boost.tar.gz boost/
mv boost.tar.gz ../../../external
cd ../..

#### cleanup
rm -rf boost_${VERSION_}
