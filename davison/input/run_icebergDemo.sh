#!/bin/bash
#$ -l h_rt=05:00:00
#$ -l rmem=4G

# add mpi modules
module load dev/PGI-compilers/17.5

# add libraries to path
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/packages/dev/pgi/linux86-64/17.5/lib"
export LDFLAGS="-rpath /usr/local/packages/dev/pgi/linux86-64/17.5/lib "$LDFLAGS

# execution bit:
./mitgcmuv > output.txt