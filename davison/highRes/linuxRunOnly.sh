#!/usr/bin/env bash
set -e

echo "Already build, clean up run directory, then make simlinks and run"
cd run
touch test.txt #ensure the directory is not empty
rm *
ln -s ../input/* .
cp ../build/mitgcmuv .
pwd
./mitgcmuv
echo "remove simlinks"
#find . -maxdepth 1 -type l -delete
echo "copy over figures and data though"
#cp ../input/*.png . 
#cp ../input/*.mat .
echo "End of script..."
