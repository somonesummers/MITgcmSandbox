#!/usr/bin/env bash
module load gcc  mvapich2
export MPI_GCC_DIR=/usr/local/pace-apps/spack/packages/linux-rhel9-x86_64_v3/gcc-12.3.0/mvapich2-2.3.7-1-qv3gjagtbx5e3rlbdy6iy2sfczryftyt
export MPI_INC_DIR=$MPI_GCC_DIR/include
export PATH=$MPI_GCC_DIR/bin:$PATH

cd build
rm *
../../../tools/genmake2 -mods ../code -mpi -rootdir ../../.. -optfile ../linux_amd64_gfortran
