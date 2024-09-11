#!/bin/bash
#SBATCH --job-name=psummers8MITgcm                     # Job name
#SBATCH --account=gts-arobel3                   # charge account
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=1G                        # Memory per core
#SBATCH --time=00:10:00                         # Duration of the job 
#SBATCH -qinferno                               # QOS Name
#SBATCH -oReport-%j.out                         # Combined output and error messages file
#SBATCH --mail-type=END,FAIL                         # Mail preferences
#SBATCH --mail-user=psummers8@gatech.edu        # E-mail address for notifications
cd $SLURM_SUBMIT_DIR                            # Change to working directory

bash linuxBuildRun.sh
