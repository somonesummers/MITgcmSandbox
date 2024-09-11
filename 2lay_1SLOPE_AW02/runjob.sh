#!/bin/bash
#SBATCH -A snic2021-1-2
#SBATCH -J 2l1Sl02
#SBATCH -t 24:00:00
#SBATCH -N 2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jonathan.wiskandt@misu.su.se
#########################################################################################

### run MITgcm and save output
rm slurm*
rm *.data
rm *.meta
rm STD*

mpprun mitgcmuv

# Script ends here
