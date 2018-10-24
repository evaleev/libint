#!/bin/sh

set -e

aclocal -I lib/autoconf
autoconf

rm -rf autom4te.cache aclocal.m4


