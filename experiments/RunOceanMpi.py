#!/usr/bin/env python3
"""
Created on Wed Mar 26 2025

@author: psummers8

Main Script for running mélange/MITgcm coupled model.
Depends on advectBergs.py, this should be incorperated. 
Depends on glaciome1d class. This should stand alone.



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
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import constants, glaciome
import glob
import pickle
import time
import copy
import signal

# resetStart = True
# freshStart = False
# icebergCoverLambda = .8
# iterationsToRun = 100
#Define if new run or not, and packing fraction of melange
sys.path.append('.')

# Import run specific settings from the experiment directory
forceMelange = False
from melangeModel import *

couplingTime = 12 * 3600

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

# This used to reset the directory to initial conditions, currently out of date, doesnt work for spun up starts
if(resetStart):
    sysPrint("Resetting this experiment directory, DELETING FILES AND FIGS. Steps to run: %i..." %iterationsToRun)
    time.sleep(1)
    # Reset iceberg files (don't need to do this so long as mélange is roughly similar)
    # os.system('rm input/icebergs_length.bin input/icebergs_depths.bin input/icebergs_widths.bin')
    # os.system('cp input/icebergs_length_init.bin input/icebergs_length.bin')
    # os.system('cp input/icebergs_depths_init.bin input/icebergs_depths.bin')
    # os.system('cp input/icebergs_widths_init.bin input/icebergs_widths.bin')

    # Clean out itermediate save states and figures
    os.system("find couplingResults/MITgcmRun_*.pickle ! -name 'MITgcmRun_00000.pickle' -type f -exec rm {} +")
    os.system("rm couplingResults/*.txt")
    os.system("rm input/pickup*.data input/pickup*.meta")
    os.system("rm figs/*")

    # Reset data file with correct start/stop times
    for line in fileinput.input('input/data'):
            if "nIter0=" in line:
                startIter = float(line[8:-2])
            elif "endTime=" in line:
                endTime = float(line[9:-2])
    sysPrint('\tadjust start iteration %i to %i' %(int(startIter),int(0)))
    replaceAll('input/data','nIter0=%i' %(int(startIter)), 'nIter0=%i' % int(0))
    # So ptracers Iter0 should actually be 0, it is when the tracer exp starts, not what to load
    # replaceAll('input/data.ptracers','Iter0=%i' %(int(startIter)), 'Iter0=%i' % int(0))
    sysPrint('\tadjust end time %i to %i' %(int(endTime),int(couplingTime)))
    replaceAll('input/data','endTime=%i' %(int(endTime)), 'endTime=%i' %(int(couplingTime)))

    #Set melange to length consistent with initial glaciome1d size
    os.system("~/.conda/envs/MITgcm/bin/python advectBergs.py >> couplingResults/out.txt")

    #Run intial MITgcm, this resets the results folder
    os.system('bash ../makeRun.sh')
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
    os.system('bash ../makeRunMpi.sh > MITgcmInitOut.txt')

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
        dataTemp = mds.rdmds("results/%s"%(prefixes[k]), endIter)
        mds.wrmds('results/%s' %prefixes[k],dataTemp,itr=endIter, dataprec='float32')
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

    if(False and index > 1 and index % 100 == 0): #block this for now
        # Reset prtracers every 100 steps, avoid saturation.
        # We read old value, then replace it with new in ptracer file.
        # We also must delete the ptracer pickup files to start at tracers=0.
        sysPrint('\t\t PTRACERS RESET every 100 iterations')
        for line in fileinput.input('input/data.ptracers'):
            if "Iter0" in line:
                oldStartIterPt = int(line[16:-1])# This assumes this is last line in data.prtacers. make -2 if more lines after
        replaceAll('input/data.ptracers','PTRACERS_Iter0=%i' %(int(oldStartIterPt)), 'PTRACERS_Iter0=%i' % int(newStartTime/dt))
        os.system('rm results/pickup_ptracers.%010i.*' %int((newStartTime)/dt))

    # We adjust the icebergs to the new mélange geometry. This could be within this script.
    os.system("~/.conda/envs/MITgcm/bin/python advectBergs.py >> couplingResults/out.txt")

    os.chdir("results")
    os.system('srun ./mitgcmuv') # srun has no outputs, all in STDOUT/STDERR.*.*
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
        dataTemp = mds.rdmds("results/%s"%(prefixes[k]), endIter)
        mds.wrmds('results/%s' %prefixes[k],dataTemp,itr=endIter, dataprec='float32')
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
