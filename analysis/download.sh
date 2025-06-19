#!/usr/bin/env bash
# set -e

doResults=true
while getopts :hr opt; do
    case $opt in 
        h) echo "help not defined"; exit ;;
        r) doResults=false ;;
        :) echo "Missing argument for option -$OPTARG"; exit 1;;
       \?) echo "Unknown option -$OPTARG"; exit 1;;
    esac
done

shift $(( OPTIND - 1 ))

if [ $# -lt 1 ]; then
  echo 1>&2 "$0: not enough arguments, need to specifiy experiment name"
  exit 2
fi

for EXP in "$@"
do
  echo $EXP
  if [ ! -d "$EXP" ]; then
    echo "$EXP does not exist locally, mkdir and downloading."
    mkdir $EXP
    mkdir $EXP/figs
    mkdir $EXP/results
    mkdir $EXP/input
  else
    echo " $EXP exists locally, syncing from PACE"
  fi
  
  #Sometimes one login is slower, can flip between here
  echo ' output file (if there)'
  rsync -azh --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/output* $EXP/
  echo ' input/'
  rsync -azh --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/input/ $EXP/input
  echo ' couplingResults/'
  rsync -azh --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/couplingResults/ $EXP/couplingResults
  echo 'figs/'
  rsync -azh --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/figs/ $EXP/figs
  if $doResults; then
    echo ' results/'
    rsync -azh --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/results/ $EXP/results
  fi
  # rsync -ah --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/output* $EXP/
  # rsync -ah --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/input/ $EXP/input
  # rsync -ah --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/results/ $EXP/results

  # echo "===== Running plotting now ====="
  # bash plotAll.sh $EXP
done
