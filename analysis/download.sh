#!/usr/bin/env bash
set -e

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
    echo "$EXP exists locally, syncing from PACE"
  fi
  
  #Sometimes one login is slower, can flip between here

  rsync -ah --info=progress2 psummers8@login-phoenix-rh7.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/output* $EXP/
  rsync -ah --info=progress2 psummers8@login-phoenix-rh7.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/input/ $EXP/input
  rsync -ah --info=progress2 psummers8@login-phoenix-rh7.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/results/ $EXP/results


  # rsync -ah --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/output* $EXP/
  # rsync -ah --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/input/ $EXP/input
  # rsync -ah --info=progress2 psummers8@login-phoenix-rh9.pace.gatech.edu:~/MITgcmSandbox/experiments/$EXP/results/ $EXP/results

  # echo "===== Running plotting now ====="
  # bash plotAll.sh $EXP
done