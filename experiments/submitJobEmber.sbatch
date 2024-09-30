#!/bin/bash
#SBATCH --job-name=Berg95                       # Job name
#SBATCH --account=gts-arobel3                   # charge account
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00                         # Duration of the job 
#SBATCH -qembers                                # QOS Name
#SBATCH -oReport-%j.out                         # Combined output and error messages file
#SBATCH --mail-type=END,FAIL                         # Mail preferences
#SBATCH --mail-user=psummers8@gatech.edu        # E-mail address for notifications
cd $SLURM_SUBMIT_DIR                            # Change to working directory

bash linuxInputOnly.sh
bash linuxBuildOnly.sh
bash linuxRunOnly.sh

#bash linuxBuildOnly.sh
