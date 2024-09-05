#!/bin/bash

# load modules (copy into command line)
module load dev/PGI-compilers/17.5

## Compile code
mkdir ./build # create build directory
cd ./build # move to build directory
# start compiling
../../tools/genmake2 -mods=../code -optfile=../../optfile/linux_ia32_pgf77_icebergmpi_netcdf_sharc -rootdir=../../
make depend
make

# Setup run directory
cd ../ # Move back to main directory
rm -r ./run
mkdir ./run
cd ./run
cp -r ../input/* . # copy contents of input to run
cp ../build/mitgcmuv . # copy MITgcm executable to run
chmod 755 ./mitgcmuv # ensure we have permissions to run MITgcm

# Run
qsub run_icebergDemo.sh # submit as batch job