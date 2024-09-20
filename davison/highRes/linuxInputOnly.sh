#!/usr/bin/env bash
set -e

module load matlab

cd input
echo 'removing old iceberg inputs'
#find . -name iceberg_* -delete
#find . *.bin -delete
#find . *.png -delete
#find . berg_generation_data.mat -delete

echo 'making new berg input'
matlab -nodisplay < gendata.m

echo "berg input done"
