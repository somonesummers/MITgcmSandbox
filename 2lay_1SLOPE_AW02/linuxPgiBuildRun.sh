#!/usr/bin/env bash
set -e

bash linuxPgiBuildOnly.sh

echo "running post build"
cd run
../build/mitgcmuv > output.txt

