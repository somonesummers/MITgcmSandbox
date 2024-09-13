#!/usr/bin/env bash
set -e
echo "Building already done, just new runtime options"
cd run
echo "clean up run folder, then make simlinks and run"
rm -r *
mkdir figs
ln -s ../input/* .
echo  'running at directory '$(pwd)
time mpirun -np 8 ../build/mitgcmuv
afplay /System/Library/Sounds/Funk.aiff
#python ../quicklook.py
#python ../visualize.py
echo -e "End of script...\n"
