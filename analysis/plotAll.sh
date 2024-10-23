#!/usr/bin/env bash
set -e

set -e

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments, need to specifiy experiment name"
  exit 2
fi

for EXP in "$@"
do
  echo $EXP
  cd $EXP

  echo '== Side View ==' 
  python ../quicklookSide.py

  echo '== End View ==' 
  python ../quicklookXSlice.py

  echo '== Map View ==' 
  python ../quicklookMap.py

  #echo '== Melt View ==' #This is now inside the cross sections as Berg Melt is 3D
  #python ../quickMeltRate.py

  cd ..
  #echo 'cleaning up pngs'
  #rm *.png
    
done
