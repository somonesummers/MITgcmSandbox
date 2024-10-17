#!/usr/bin/env bash
set -e

set -e

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments, need to specifiy experiment name"
  exit 2
fi

for EXP in "$@"
do
  cd $EXP

  echo '== Side View ==' 
  python ../quicklookSide.py

  echo '== End View ==' 
  python ../quicklookXSlice.py

  echo '== Map View ==' 
  python ../quicklookMap.py

  echo '== Melt View ==' 
  python ../quickMeltRate.py

  #echo 'cleaning up pngs'
  #rm *.png
    
done
