#!/bin/sh

set -e

sudo add-apt-repository ppa:boost-latest/ppa -y
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y 
sudo apt-get update
sudo apt-get install -y g++-$GCC_VERSION llvm-3.4 llvm-3.4-dev libgmp-dev libeigen3-dev libboost1.55-dev
