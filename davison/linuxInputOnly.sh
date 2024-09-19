#!/usr/bin/env bash
set -e

cd input
echo 'Make new input files'
matlab -nodisplay < ../matlabFunctions/gendata.m
echo 'Done making input'
