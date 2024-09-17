#!/usr/bin/env bash
export PKG_CONFIG_PATH=/usr/local/pace-apps/spack/packages/linux-rhel9-x86_64_v3/gcc-12.3.0/mvapich2-2.3.7-1-qv3gjagtbx5e3rlbdy6iy2sfczryftyt/bin/mpicc
cd build
echo -e "\nCleaning up previous builds...\n"
make Clean
../../tools/genmake2 -mods ../code -mpi -optfile ../../tools/build_options/linux_amd64_gfortran
echo -e "\nDone compiling, moving to make depend...\n"
make depend -s
echo -e "\nDone with make depend, moving to make...\n"
make -s
cd ../run
echo -e "\nclean up run folder, then make simlinks and run"
rm -r *
mkdir figs
ln -s ../input/* .
echo  'running at directory '$(pwd)
time srun ../build/mitgcmuv
python ../visualize.py
echo "End of script..."
