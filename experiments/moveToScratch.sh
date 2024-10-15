#!/usr/bin/env bash
set -e

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments"
  exit 2
fi

for NAME in "$@"
do
	mkdir ../../scratch/exps/$NAME
	rsync -ah --info=progress2 $NAME ../../scratch/exps/
	rm -r $NAME
	ln -s ../../scratch/exps/$NAME $NAME
done
