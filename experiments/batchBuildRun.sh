#!/usr/bin/env bash
set -e

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments, specify experiments to run"
  exit 2
fi

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     MACHINE="Linux";;
    Darwin*)    MACHINE="Mac";;
esac
echo "Idenitfied machine as ${MACHINE}"


for NAME in "$@"
do
	echo "move in directory $NAME..."
	cd $NAME
	if [ "$MACHINE" == "Linux" ];
	then
		bash ../makeBuild.sh ../../.. -mpi
		echo "submitting job"
		sbatch submitBatch.sh
		echo "moving back"
	else
		bash ../makeBuild.sh ../../..
		bash ../makeRun.sh
		echo "moving back"
	fi
	cd ..
done

if [ "$MACHINE" == "Mac" ];
then
   afplay /System/Library/Sounds/Funk.aiff &
fi 

echo "(⌐■_■) Done with batch"
