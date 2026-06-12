#!/usr/bin/env python3
"""
Created on Wed Mar 26 2025

@author: psummers8

Main Script for running iceberg advection run




"""

import numpy as np
from matplotlib import pyplot as plt
import sys
import fileinput
import shutil
from MITgcmutils import mds
import os
# sys.path.append('/hdd/glaciome/models/glaciome1D')
# sys.path.insert(0, '')
import time
import copy
import signal
import platform
import re

# resetStart = True
# freshStart = False
# icebergCoverLambda = .8
# iterationsToRun = 100
#Define if new run or not, and packing fraction of melange
sys.path.append('.')

# Import run specific settings from the experiment directory
forceMelange = False
from melangeModel import *

couplingTime = 6 * 3600

# Method to return to main output as well as status text file
def sysPrint(stringIn): 
        print(stringIn)
        os.system('echo ""' + stringIn + '"" >> couplingResults/out.txt')

# Edits files in place, used for adjusting MITgcm files
def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

# Track exceptions
mitgcmWarningFlag = 1 #if it fails first run, should error
mitgcmWarningCount = 0

OSX = platform.system()
current_directory = os.getcwd()
baseDir = ''
CondaDir = ''
if OSX == 'Darwin':
    CondaDir = '~/anaconda3/envs/MITgcm/bin/'
elif "cluster" in current_directory:
    CondaDir = '~/.conda/envs/MITgcm/bin/'
    OSX = 'Tufts'
elif ("storage" in current_directory) or ("hcoda" in current_directory):
    CondaDir = '~/.conda/envs/MITgcm/bin/'
    OSX = 'PACE'
else:
    print(current_directory)
    raise Exception("unknown OSX or running location, please configure")

sysPrint(f'Operating system identified as {OSX}')

nPx = 0 
nPy = 0
isMPI = False
for line in fileinput.input('code/SIZE.h'):
    if "nPx =" in line:
        nPx = float(re.sub(r"\D", "", line)) #strips all non digits out of line
    if "nPy =" in line:
        nPy = float(re.sub(r"\D", "", line))  #strips all non digits out of line
if nPx * nPy > 1:
    isMPI = True
if(nPy*nPx == 0):
    raise Exception("Cannot identify cores required, please configure")
else:
    sysPrint(f'Identified cores as {nPx} x {nPy} = {nPy*nPx}, isMPI = {isMPI}')



# This used to reset the directory to initial conditions, currently out of date, doesnt work for spun up starts
if(resetStart):
    sysPrint("Not Supported")

elif(freshStart): #This distinguises between a new coupled run, or continuing a paused coupled run
    #Make couplingResults directory if not there already
    os.system("mkdir -p couplingResults")
    os.system("rm -f couplingResults/out.txt")
    sysPrint('Running from fresh start...')
    sysPrint("Assuming fresh MITgcm directory. Steps to run: %i..." %iterationsToRun)
    time.sleep(1)
    #prime the iceberg files
    #os.system('cp input/icebergs_length_init.bin input/icebergs_length.bin')
    #os.system('cp input/icebergs_depths_init.bin input/icebergs_depths.bin')
    #os.system('cp input/icebergs_widths_init.bin input/icebergs_widths.bin')
    
    # run initial MITgcm
    #mpirun flags, these help reduce minor errors and unneeded output
    mpi_flags = ''
    if(OSX == 'Tufts'):
        mpi_flags = '--quiet --mca btl_vader_single_copy_mechanism none'
    elif(OSX == 'PACE'):
        mpi_flags = ''
    elif(OSX == 'Darwin'):
        mpi_flags = ''
    # run initial MITgcm
    os.system("touch results/test.txt")
    os.chdir("results")
    os.system("rm *")
    os.system("cp ../build/mitgcmuv .")
    os.system("ln -s ../input/* .")
    if(isMPI):
        if(OSX == 'Tufts'):
            os.system(f"mpirun {mpi_flags} -n {nPy*nPx} ./mitgcmuv")
        elif(OSX == 'PACE'): 
            os.system(f"srun {mpi_flags} ./mitgcmuv") #srun is smart and allocated cores automatically
        elif(OSX == 'Darwin'):
            os.system(f"mpirun {mpi_flags} -n {nPy*nPx} ./mitgcmuv")
        else:
            raise ExceptionType("unknown OSX or running location, please configure")
    else:
        os.system(f"./mitgcmuv > MITgcm_output.txt")
    os.chdir("../")

    # Clean tile level files
    dt = 0.0
    oldStartIter = 0   
    for line in fileinput.input('input/data'):
            if "deltaT=" in line:
                dt = float(line[8:-2])
            if "nIter0" in line:
                oldStartIter = int(line[8:-2])

    prefixes = ['BRGFlx','dynDiag','ptraceDiag','plumeDiag']
    endIter = int((couplingTime)/dt)
    sysPrint('\tcondensing diagnostic tile files to global files iter:%i' %endIter)
    for k in range(len(prefixes)):
        sysPrint('\t\t == %s ==' %prefixes[k])
        dataTemp, it, meta = mds.rdmds("results/%s"%(prefixes[k]), endIter, returnmeta=True)
        mds.wrmds('results/%s' %prefixes[k],dataTemp,itr=endIter, dataprec='float32', fields=meta['fldlist'])
        os.system('rm results/%s.%010i.0*.0*' %(prefixes[k],endIter)) #picks out tile level files


else:
    # This means were starting up an existing coupled run. Surprisingly easy for us to
    # Changing MITgcm dt is possible here, but breaks lots of later plotting so dont do it lightly
    sysPrint('Running from existing states...')
    time.sleep(1)


for ii in range(iterationsToRun):
    # Welcome to the loop. We assume there is a MITgcm file at this point, and a mélange geometry existing.
    # ii is NOT the final iteration to run, simply just counting up. ii Never actually used. 
    # iteration index comes from glaciome pickle file names and MITgcm time index.

    start_time = time.time() #time steps, nice to know

    #Find time step to take from MITgcm files
    maxStep = 0
    for file in os.listdir('results'):
        # print(file)
        if "BRGFlx.0" in file:
            words = file.split(".")
            # print(words[1])  
            if int(words[1]) > maxStep:
                maxStep = int(words[1])

    # As mentioned, this script can handle chaning dt, but plotting cant currently
    dt = 0.0
    oldStartIter = 0   
    for line in fileinput.input('input/data'):
            if "deltaT=" in line:
                dt = float(line[8:-2])
            if "nIter0" in line:
                oldStartIter = int(line[8:-2])
    # print('\tdt is loaded as', dt)

    #Load all our MITgcm data
    x = np.squeeze(mds.rdmds("results/XC")[0,:])
    y = np.squeeze(mds.rdmds("results/YC")[:,0])
    z = np.squeeze(mds.rdmds("results/RC"))
    dx = 2*x[0]
    dy = 2*y[0]
    ny = len(y)
    nx = len(x)
    nz = np.squeeze(np.shape(z)[0])

    dataMITgcm = mds.rdmds("results/BRGFlx", maxStep)
    openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((nz,ny,nx))

    index = int(maxStep)*dt/(couplingTime)

    sysPrint('\tPrepping MITgcm for step %.2f' %(index+1))
    newStartTime = (index)*couplingTime
    oldStartTime = (oldStartIter)*dt #read directly from input/data now
    # we delete the old pickups to reduce file size and file count. We could strategically leave them every
    # 100 iterations or so, but currently we do not. 
    os.system('rm results/pick*.%010i.*' %int((oldStartTime-couplingTime)/dt)) #save last one, but delete 2 ago

    sysPrint('\tadjust start iteration %i to %i' %(oldStartIter,int(newStartTime/dt)))
    replaceAll('input/data','nIter0=%i' %(int(oldStartIter)), 'nIter0=%i' % int(newStartTime/dt))
    sysPrint('\tadjust end time %i to %i' %(int(newStartTime),int(newStartTime+couplingTime)))
    replaceAll('input/data','endTime=%i' %(int(newStartTime)), 'endTime=%i' %(int(newStartTime+couplingTime)))

    # Advect icebergs
    os.system(f"{CondaDir}python advectBergsOcean.py >> couplingResults/out.txt")

    if(OSX == 'Tufts'):
        mpi_flags = '--quiet --mca btl_vader_single_copy_mechanism none'
    elif(OSX == 'PACE'):
        mpi_flags = ''
    elif(OSX == 'Darwin'):
        mpi_flags = ''
    # run on existing MITgcm files, no clearing of the folder
    os.chdir("results")
    if(isMPI):
        if(OSX == 'Tufts'):
            os.system(f"mpirun {mpi_flags} -n {nPy*nPx} ./mitgcmuv")
        elif(OSX == 'PACE'): 
            os.system(f"srun {mpi_flags} ./mitgcmuv") #srun is smart and allocated cores automatically
        elif(OSX == 'Darwin'):
            os.system(f"mpirun {mpi_flags} -n {nPy*nPx} ./mitgcmuv")
        else:
            raise ExceptionType("unknown OSX or running location, please configure")
    else:
        os.system(f"./mitgcmuv > ../MITgcm_output.txt")
    os.chdir("../")
    
    # Create global files to reduce file counts
    # MITgcm MPI saves outputs per tile, per timestep. For 20 cores this makes 40 files per diagnostic
    # per timestep, which quickly become insane for long runs. This step collects all those tile files into
    # 2 global files (*.data, *.meta), then deletes the tiles files. This reduces files count by a factor of 20
    # This VASTLY improves data transfer and compression speeds. 
    prefixes = ['BRGFlx','dynDiag','ptraceDiag','plumeDiag']

    endIter = int((newStartTime + couplingTime)/dt)
    sysPrint('\tcondensing diagnostic tiles files to global files iter:%i' %endIter)
    for k in range(len(prefixes)):
        sysPrint('\t\t == %s ==' %prefixes[k])
        dataTemp, it, meta = mds.rdmds("results/%s"%(prefixes[k]), endIter, returnmeta=True)
        mds.wrmds('results/%s' %prefixes[k],dataTemp,itr=endIter, dataprec='float32', fields=meta['fldlist'])
        os.system('rm results/%s.%010i.0*.0*' %(prefixes[k],endIter))
    
    sysPrint('\t\tSeconds to run coupled step: %.4f' % (time.time() - start_time))

    # Now start all over again, wohooo I love loops

## After final iteration, we can clean up the tile files to global files, mostly for conistancy
prefixes = ['Depth','DXC','DXF','DXG','DXV','DYC','DYF','DYG','DYU','hFacC','hFacS','hFacW',
                    'maskInC','maskInS','maskInW','RAC','RAS','RAW','RAZ','XC','XG','YC','YG']
sysPrint('\tcondensing grid tile files to global files')
for k in range(len(prefixes)):
    sysPrint('\t\t== %s ==' %prefixes[k])
    dataTemp = mds.rdmds("results/%s"%(prefixes[k]))
    mds.wrmds('results/%s' %prefixes[k],dataTemp,dataprec='float32')
    os.system('rm results/%s.0*.0*' %(prefixes[k])) #picks out tile level files
