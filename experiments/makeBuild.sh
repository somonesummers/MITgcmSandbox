#!/usr/bin/env bash
set -e

#module load nvhpc/24.5
#export PGI=/usr/local/pace-apps/manual/packages/nvhpc/24.5
#export PATH=$PGI/linux86-64/24.5/compilers/bin:$PATH
#export MANPATH=$MANPATH:$PGI/linux86-64/24.5/compilers/man
#export LM_LICENSE_FILE=$PGI/license/LICENSE.txt
if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments"
  exit 2
fi

ROOT=$1

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     MACHINE="Linux";;
    Darwin*)    MACHINE="Mac";;
esac
echo "Idenitfied machine as ${MACHINE}"

cd build
echo " Cleaning up previous builds..."
make Clean
if [ "$MACHINE" == "Linux" ];
then
	#BUILD_FILE='linux_amd64_pgf77_pace'
	BUILD_FILE='linux_amd64_gfortran'
else
	BUILD_FILE='darwin_amd64_gfortran'
fi

echo " Building with ${BUILD_FILE} at " $(pwd)
$ROOT/tools/genmake2 -mods ../code -optfile $ROOT/tools/build_options/$BUILD_FILE -rootdir $ROOT
echo " Done compiling, moving to make depend..."
make depend -s
echo " Done with make depend, moving to make..."
make -s
echo "Done building..."
