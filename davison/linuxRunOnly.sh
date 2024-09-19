#!/usr/bin/env bash
set -e

echo "Already build, clean up run folder, then make simlinks and run"
cd run
touch test.txt #this ensures the dir is not empty
rm *
ln -s ../input/* .
cp ../build/mitgcmuv .
pwd
./mitgcmuv > output.txt
echo "remove simlinks"
find /path/to/directory -maxdepth 1 -type l -delete
echo "copy over figures though"
cp ../input/*.png . 
cp ../input/*.mat .
echo "end of script..."
