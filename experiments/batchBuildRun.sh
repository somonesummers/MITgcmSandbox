#!/usr/bin/env bash
set -e

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments, specify experiments to run"
  exit 2
fi

for NAME in "$@"
do
	echo "move in directory $NAME..."
	cd $NAME
	bash ../makeBuild.sh ../../.. -mpi
	echo "submitting job"
	sbatch submitBatch.sh
	echo "moving back"
	cd ..
done
