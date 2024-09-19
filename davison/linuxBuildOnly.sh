#!/usr/bin/env bash
set -e

#module load nvhpc/24.5
#export PGI=/usr/local/pace-apps/manual/packages/nvhpc/24.5
#export PATH=$PGI/linux86-64/24.5/compilers/bin:$PATH
#export MANPATH=$MANPATH:$PGI/linux86-64/24.5/compilers/man
#export LM_LICENSE_FILE=$PGI/license/LICENSE.txt

cd build
echo " Cleaning up previous builds..."
make Clean
#BUILD_FILE='linux_amd64_pgf77_pace'
BUILD_FILE='linux_amd64_gfortran'
echo " Building with ${BUILD_FILE}"
../../tools/genmake2 -mods ../code -optfile ../../tools/build_options/$BUILD_FILE -rootdir ../..
echo " Done compiling, moving to make depend..."
make depend -s
echo " Done with make depend, moving to make..."
make -s
cd ../run
echo " Clean up run folder, then make simlinks and run"
rm *
ln -s ../input/* .
cp ../build/mitgcmuv .
echo "End of script, run directory is configured..."
