#!/bin/bash
#SBATCH --job-name=shelf500eq                     # Job name
#SBATCH --account=gts-arobel3                   # charge account
#SBATCH -N1 --ntasks-per-node=5
#SBATCH --mem=2G
#SBATCH --time=01:00:00                         # Duration of the job 
#SBATCH -qinferno                                # QOS Name (embers is free, inferno is paid)
#SBATCH -oReport-%j.out                         # Combined output and error messages file
#SBATCH --mail-type=END,FAIL                         # Mail preferences
#SBATCH --mail-user=psummers8@gatech.edu        # E-mail address for notifications
cd $SLURM_SUBMIT_DIR                            # Change to working directory
set -e
#bash ../makeInput.sh
bash ../makeBuild.sh ../../.. -mpi
bash ../makeRunMpi.sh
