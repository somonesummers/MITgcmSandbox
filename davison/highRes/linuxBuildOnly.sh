#!/usr/bin/env bash
set -e

cd build
echo " Cleaning up previous builds..."
#make Clean
#BUILD_FILE='linux_amd64_pgf77_pace'
BUILD_FILE='linux_amd64_gfortran'
echo " Building with ${BUILD_FILE}"
../../../tools/genmake2 -mods ../code -optfile ../../../tools/build_options/$BUILD_FILE -rootdir ../../..
echo " Done compiling, moving to make depend..."
make depend -s
echo " Done with make depend, moving to make..."
make -s
echo "Building done"
