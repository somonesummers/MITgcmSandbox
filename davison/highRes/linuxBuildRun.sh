#!/usr/bin/env bash
set -e
bash linuxBuildOnly.sh
cd run
echo "clean up run folder, then make simlinks and run"
rm *
ln -s ../input/* .
ln -s ../build/mitgcmuv .
pwd
./mitgcmuv > output.txt
echo "End of script..."
