#!/bin/bash
#SBATCH --job-name=CHANGEME                     # Job name
#SBATCH --account=gts-arobel3                   # charge account
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=2G   
#SBATCH --time=08:00:00                         # Duration of the job 
#SBATCH -qembers/inferno                                # QOS Name (embers is free, inferno is paid)
#SBATCH -oReport-%j.out                         # Combined output and error messages file
#SBATCH --mail-type=END,FAIL                         # Mail preferences
#SBATCH --mail-user=psummers8@gatech.edu        # E-mail address for notifications
cd $SLURM_SUBMIT_DIR                            # Change to working directory
set -e

#bash ../makeBuild.sh ../../..
#bash ../makeRun.sh

bash ../makeBuild.sh ../../.. -mpi
bash ../makeRunMpi.sh
