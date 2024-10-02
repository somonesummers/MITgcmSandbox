#!/bin/bash
#SBATCH --job-name=B95_250_17280               # Job name
#SBATCH --account=gts-arobel3                   # charge account
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00                         # Duration of the job 
#SBATCH -qembers                                # QOS Name
#SBATCH -oReport-%j.out                         # Combined output and error messages file
#SBATCH --mail-type=END,FAIL                         # Mail preferences
#SBATCH --mail-user=psummers8@gatech.edu        # E-mail address for notifications
cd $SLURM_SUBMIT_DIR                            # Change to working directory
set -e
bash ../makeInput.sh
bash ../makeBuild.sh ../../..
bash ../makeRun.sh
