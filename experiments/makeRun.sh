#!/usr/bin/env bash
set -e

echo "Already build, clean up run folder, then make simlinks and run"
cd run
touch test.txt #this ensures the dir is not empty
rm *
ln -s ../input/* .
cp ../build/mitgcmuv .
./mitgcmuv > output.txt

echo "Done running..."
