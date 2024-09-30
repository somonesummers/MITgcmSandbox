#!/usr/bin/env bash
set -e

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     MACHINE="Linux";;
    Darwin*)    MACHINE="Mac";;
esac
echo "Identify machine as ${MACHINE}"
	
cd input
echo 'Make new input files'

if [ "$MACHINE" == "Linux" ];
then
	module load matlab
	matlab -nodisplay < gendata.m
else
	/Applications/MATLAB_R2024a.app/bin/matlab -nodisplay < gendata.m
fi

echo 'Done making input...'
