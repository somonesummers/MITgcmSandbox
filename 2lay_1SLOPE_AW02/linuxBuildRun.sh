#!/usr/bin/env bash
modeule load gcc mvapich2
cd build
echo " Cleaning up previous builds..."
make Clean
/../../tools/genmake2 -mods ../code -optfile ../../tools/build_options/linux_amd64_gfortran
echo "Done compiling, moving to make depend..."
make depend -s
echo "Done with make depend, moving to make..."
make -s
cd ../run
echo "clean up run folder, then make simlinks and run"
rm -r *
mkdir figs
ln -s ../input/* .
rm mitgcmuv
ln -s ../build/mitgcmuv .
pwd
./mitgcmuv > output.txt
echo "End of script..."
