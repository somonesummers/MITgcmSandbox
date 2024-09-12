#!/usr/bin/env bash
# export MPI_HOME=/abd/whereIam #if it can't find MPI itself...
cd build
echo -e "\nCleaning up previous builds...\n"
make Clean
../../tools/genmake2 -mods ../code -mpi -optfile ../../tools/build_options/linux_amd64_gfortran
echo -e "\nDone compiling, moving to make depend...\n"
make depend -s
echo -e "\nDone with make depend, moving to make...\n"
make -s
cd ../run
echo -e "\nclean up run folder, then make simlinks and run"
rm -r *
mkdir figs
ln -s ../input/* .
echo  'running at directory '$(pwd)
time srun ../build/mitgcmuv
afplay /System/Library/Sounds/Funk.aiff
python ../visualize.py
echo "End of script..."