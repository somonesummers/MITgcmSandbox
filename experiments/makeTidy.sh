#!/usr/bin/env bash
set -e

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments"
  exit 2
fi

RUN_NAME=$1
mkdir ../../analysis/$RUN_NAME
mkdir ../../analysis/$RUN_NAME/input
FINAL_DIR=../../analysis/$RUN_NAME
echo $FINAL_DIR
cd run

echo "Cleaning up run folder, remove simlinks"
find . -maxdepth 1 -type l -delete
echo "Copy over non-simlinks"
cp * ../$FINAL_DIR
echo "Copy over input data"
cp ../input/* ../$FINAL_DIR/input
echo "Clear run folder"
rm *
echo "Done tidying"