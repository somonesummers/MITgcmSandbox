#!/usr/bin/env bash
set -e

TIME="$(date +"%H%M%S")" 
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     MACHINE="Linux";;
    Darwin*)    MACHINE="Mac";;
esac
echo "Idenitfied machine as ${MACHINE}"
echo "Already build, clean up run folder, then make simlinks and run"
cd results
touch test.txt #this ensures the dir is not empty
rm *
ln -s ../input/* .
cp ../build/mitgcmuv .

if [ "$MACHINE" == "Mac" ];
then
	if [ $# -lt 1 ]; then
  		echo 1>&2 "$0: not enough arguments, need to specify number of cores"
  		exit 2
	fi
	mpirun -np 5 ./mitgcmuv > Report$TIME.txt
else
	srun ./mitgcmuv
fi

echo "Done running..."
