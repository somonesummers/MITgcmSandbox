#!/bin/bash
#$  -m eba   -M bd41@st-andrews.ac.uk
#$ -l h_rt=50:00:00
#$ -l rmem=4G
#$ -pe mpi 20
#
# add mpi modules
module load dev/PGI-compilers/17.5
module load mpi/openmpi/2.0.1/pgi-17.5

# add libraries to path
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/packages/dev/pgi/linux86-64/17.5/lib"
export LDFLAGS="-rpath /usr/local/packages/dev/pgi/linux86-64/17.5/lib "$LDFLAGS
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/packages/mpi/openmpi/2.0.1/pgi-17.5/lib"
export LDFLAGS="-rpath /usr/local/packages/mpi/openmpi/2.0.1/pgi-17.5/lib "$LDFLAGS
#
# execution bit:
# simple
#./mitgcmuv > output.txt
# MPI
mpirun ./mitgcmuv > output.txt 
# add prefix: may help to locate libraries on run nodes?
# mpirun --prefix /usr/local/mpi/pgi/openmpi/1.6.4 ./mitgcmuv > output.txt