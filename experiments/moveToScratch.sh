#!/usr/bin/env bash
set -e

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments"
  exit 2
fi

for NAME in "$@"
do
	echo "moving $NAME..."
	mkdir ../../scratch/exps/$NAME
#	cd $NAME/build
#	make Clean
#	cd ../..
	rsync -ah --info=progress2 $NAME ../../scratch/exps/
	rm -r $NAME
	ln -s ../../scratch/exps/$NAME .
done
