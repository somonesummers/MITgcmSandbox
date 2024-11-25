#!/usr/bin/env bash
# No set -e as we want every script to run

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments, need to specifiy experiment name"
  exit 2
fi

for EXP in "$@"
do
  echo $EXP
  cd $EXP

  echo '== Cross View ==' 
  python ../crossPlot.py

  echo '== Side View ==' 
  python ../quicklookSide.py

  echo '== End View ==' 
  python ../quicklookXSlice.py

  echo '== Map View ==' 
  python ../quicklookMap.py

  echo '== Depth plots ==' 
  python ../depthPlotXslice.py

  echo '== TS plot =='
  python ../quickTSplot.py

  # echo '== Plume plot =='
  # python ../quickPlumePlot.py
  
  cd ..    
done
