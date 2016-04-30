#!/bin/sh

set -e

sudo add-apt-repository 'deb http://llvm.org/apt/precise/ llvm-toolchain-precise-3.7 main'
sudo wget -O - http://llvm.org/apt/llvm-snapshot.gpg.key|sudo apt-key add -

sudo add-apt-repository ppa:boost-latest/ppa -y
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y 
sudo apt-get update

sudo apt-get install -y g++-$GCC_VERSION clang-3.7 libgmp-dev libeigen3-dev libboost1.55-dev
