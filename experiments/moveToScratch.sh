#!/usr/bin/env bash
set -e

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments"
  exit 2
fi

NAME=$1

mkdir ../../scratch/exps/$NAME
mv $NAME ../../scratch/exps/
ln -s ../../scratch/exps/$NAME $NAME
