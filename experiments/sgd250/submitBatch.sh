#!/bin/bash 
#SBATCH -J sgd250 # job name 
#SBATCH -o output_%j.txt # output and error file name (%j expands to jobID)
#SBATCH --account=gts-arobel3    #charge account
#SBATCH -N1 --ntasks-per-node=10   #total number of nodes,CPUs requested
#SBATCH --mem-per-cpu=1G
#SBATCH -qinferno
#SBATCH -t 01:08:25 # run time (hh:mm:ss)
#SBATCH --mail-user=psummers8@gatech.edu
#SBATCH --mail-type=end,fail  # email me when the job finishes/fails

cd $SLURM_SUBMIT_DIR    # Change to working directory
set -e
bash ../makeRunMpi.sh