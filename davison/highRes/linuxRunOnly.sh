#!/usr/bin/env bash
set -e

echo "Already build, clean up run folder, then make simlinks and run"
cd run
rm *
ln -s ../input/* .
cp ../build/mitgcmuv .
pwd
./mitgcmuv > output.txt
echo "End of script..."
