#!/usr/bin/env bash
cd build
echo " Cleaning up previous builds..."
make Clean
../../../../tools/genmake2 -mods ../code -optfile ../../../../tools/build_options/darwin_amd64_gfortran
echo "Done compiling, moving to make depend..."
make depend -s
echo "Done with make depend, moving to make..."
make -s
cd ../run
echo "clean up run folder, then make simlinks and run"
rm -r *
mkdir figs
ln -s ../input/* .
ln -s ../build/mitgcmuv .
echo  'running at directory '$(pwd)
./mitgcmuv > output.txt
python ../visualize.py
echo "End of script..."
