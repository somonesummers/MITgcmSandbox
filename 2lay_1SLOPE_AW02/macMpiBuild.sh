#!/usr/bin/env bash
set -e
export MPI_HOME=/opt/homebrew/
cd build
echo -e "\nCleaning up previous builds...\n"
make Clean
../../tools/genmake2 -mods ../code -mpi -optfile ../../tools/build_options/darwin_amd64_gfortran
echo -e "\nDone compiling, moving to make depend...\n"
make depend -s
echo -e "\nDone with make depend, moving to make...\n"
make -s
cd ../run
echo -e "\nclean up run folder, then make simlinks and run"
rm -r *
mkdir figs
ln -s ../input/* .
echo "End of build script"
