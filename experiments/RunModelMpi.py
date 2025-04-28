#!/usr/bin/env python3
"""
Created on Wed Mar 26 2025

@author: psummers8
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
# resetStart = True
# freshStart = False
# icebergCoverLambda = .8
# iterationsToRun = 100
#Define if new run or not, and packing fraction of melange
sys.path.append('.')

forceMelange = False
from melangeModel import *

def sysPrint(stringIn):
        print(stringIn)
        os.system('echo ""' + stringIn + '"" >> couplingResults/out.txt')

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

# This is only needed for a true fresh start
# print(resetStart,freshStart)
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
    sysPrint('\tadjust end time %i to %i' %(int(endTime),int(86400)))
    replaceAll('input/data','endTime=%i' %(int(endTime)), 'endTime=%i' %(int(86400)))

    #Set melange to length consistent with initial glaciome1d size
    os.system("~/.conda/envs/MITgcm/bin/python advectBergs.py >> couplingResults/out.txt")

    #Run intial MITgcm, this resets the results folder
    os.system('bash ../makeRun.sh')
elif(freshStart):
    #Make couplingResults directory if not there already
    os.system("mkdir -p couplingResults")
    sysPrint('Running from fresh start...')
    sysPrint("Resetting the MITgcm directory. Steps to run: %i..." %iterationsToRun)
    time.sleep(1)
    #prime the iceberg files
    #os.system('cp input/icebergs_length_init.bin input/icebergs_length.bin')
    #os.system('cp input/icebergs_depths_init.bin input/icebergs_depths.bin')
    #os.system('cp input/icebergs_widths_init.bin input/icebergs_widths.bin')
    
    #run initial MITgcm
    os.system('bash ../makeRunMpi.sh > MITgcmInitOut.txt')
else:
    sysPrint('Running from existing states...')
    time.sleep(1)


for ii in range(iterationsToRun):
    # Welcome to the loop. We assume there is a MITgcm file, and a mélange geometry existing.
    # ii is NOT the final iteration to run, simply just counting up. ii Never actually used. 

    start_time = time.time()

    #Find time step to take from MITgcm files
    maxStep = 0
    for file in os.listdir('results'):
        # print(file)
        if "BRGFlx.0" in file:
            words = file.split(".")
            # print(words[1])  
            if int(words[1]) > maxStep:
                maxStep = int(words[1])

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

    # Load mélange geometry

    files = sorted(glob.glob('couplingResults/MITgcmRun_[0-9][0-9][0-9][0-9][0-9].pickle'))
    # print(files)
    with open(files[-1], 'rb') as file:
        data = pickle.load(file)
        file.close()
        index = int(str(file).split(".")[1].split("_")[1])
        
    sysPrint('=== glaciome1D Running, load MITgcm data from end of day %.2f, using melange index %05i === ' %(int(maxStep)*dt/(24*3600),index))
    
    # make sure nothing funny is going on. Depending on when last cycle failed you may have to delete the highest numbered .pickle file
    if(int(maxStep)*dt/(24*3600) - index != 1):
        raise Exception('MITgcm (BRGFlx_) and glaciome1D (MITgcmRun_)are out of sync, glaciome index should be 1 behind MITgcm day')

    # Calculate melt rate. This is the net Freshwaterflux (sum verically, avg across)/ surface area of melange
    # It is very imporant to extend meltrate beyond mélange with value of last cell. Otherwise mélange will grow
    # continuiously. 

    lambdaHelper = 1-openFrac[0,1:-1,1:]
    lambdaHelper[lambdaHelper == 0] = 0.01
    # data is of list ['BRGfwFlx','BRGhtFlx','BRGmltRt','BRG_TauX','BRG_TauY']
    b_mitgcm = -1*np.nanmean(np.nansum(dataMITgcm[0,:,1:-1,1:],axis=0)/lambdaHelper,axis=0)/(dx*dy)*(24*3600)
    b_mitgcm[b_mitgcm == 0] = b_mitgcm[b_mitgcm != 0][-1] # We fill in the tail of 0s with the value of the melange toe

    x_mitgcm = x[1:]-x[1] #ensure starts at 0, MITgcm has a glacier for first cell

    meltHelper = b_mitgcm.copy()
    meltHelper[meltHelper==0] = np.min(meltHelper[meltHelper != 0])
    sysPrint('\tMelt Rates min/mean/max: %.2f, %.2f, %.2f [m/day]:' %(np.min(meltHelper),np.mean(meltHelper),np.max(meltHelper)))
    # plt.figure()
    # plt.plot(x_mitgcm,b_mitgcm)
    # plt.show()
    # plt.close()
    data.dt = 1.0/365 # 1 day
    data.X_externalGrid = x_mitgcm
    data.B_externalGrid = b_mitgcm*365 #days to years

    # Ready to now run GLACIOME1D. We run 1 timestep with dt set to the coupling timestep. If this creates CLF issues
    # we can run more steps with a reduced dt, or couple with MITgcm more frequently, but this isn't likely to be needed. 
    oldL = data.L
    oldH = data.H0
    if(forceMelange):
        data.Uc = 6000 + 1000 * np.sin(2 * np.pi * index/365)
        data.Ut = 6000 + 1000 * np.sin(2 * np.pi * index/365) # Must change together unless terminus moving
        sysPrint('\tNew Uc: %.3f m/year'%(data.Uc))
    sysPrint('\tPrevious Length: %.3f, H0: %.3f'%(oldL, oldH))
    # data.steadystate()
    data.prognostic(method='lm')
    sysPrint('\tNew      Length: %.3f, H0: %.3f, ∆L: %.3f, ∆H0: %.3f'%(data.L,data.H0,data.L-oldL,data.H0-oldH,))
    sysPrint('\t\tSeconds to run GLACIOME step: %.4f' % (time.time() - start_time))
    
    H = np.concatenate(([data.H0], data.H, [data.HL]))
    X = data.X
    X_ = np.concatenate(([data.X[0]],data.X_,[data.X[-1]]))
    np.save('input/melangeH',H)
    np.save('input/melangeX',X_)

    # save itermediate steps for plotting
    data.save('couplingResults/MITgcmRun_%05i.pickle' %(index+1))
    sysPrint('\tSaved as MITgcmRun_%05i.pickle' %(index+1))

    # Preparing MITgcm for next run. 
    # We adust the data file, copy the last pickup file to input directory 

    sysPrint('\tPrepping MITgcm for day %.2f' %(index+2))
    newStartTime = (index+1)*24*3600
    oldStartTime = (oldStartIter)*dt #read directly now
    sysPrint('\tmoving pickup')
    sysPrint('\tremoving old pickups, keeping last 2')
    os.system('rm results/pick*.%010i.*' %int((oldStartTime-24*3600)/dt)) #save last one, but delete 2 ago

    sysPrint('\tadjust start iteration %i to %i' %(oldStartIter,int(newStartTime/dt)))
    replaceAll('input/data','nIter0=%i' %(int(oldStartIter)), 'nIter0=%i' % int(newStartTime/dt))
    sysPrint('\tadjust end time %i to %i' %(int(newStartTime),int(newStartTime+24*3600)))
    replaceAll('input/data','endTime=%i' %(int(newStartTime)), 'endTime=%i' %(int(newStartTime+24*3600)))

    # We adjust the icebergs to the new mélange geometry. This could be within this script.

    os.system("~/.conda/envs/MITgcm/bin/python advectBergs.py >> couplingResults/out.txt")

    os.chdir("results")
    os.system('srun ./mitgcmuv') # srun has no outputs, all in STDOUT/STDERR.*.*

    # Create global files to reduce file counts
    prefixes = ['BRGFlx','dynDiag','ptraceDiag','plumeDiag']

    endIter = int((newStartTime + 24*3600)/dt)
    sysPrint('\tcondensing grid files to global files')
    for k in range(len(prefixes)):
        sysPrint('\t\t %s' %prefixes[k])
        dataTemp = mds.rdmds("%s"%(prefixes[k]), endIter)
        mds.wrmds('%s' %prefixes[k],dataTemp,itr=endIter, dataprec='float32')
        os.system('rm %s.0*.*' %(prefixes[k],endIter))
    os.chdir("../")



    dataMITgcmOcean = mds.rdmds("results/dynDiag", int((newStartTime + 24*3600)/dt))
    maxU = np.max((dataMITgcmOcean[2,:,:,:]**2 + dataMITgcmOcean[3,:,:,:]**2 + dataMITgcmOcean[2,:,:,:]**2)**(.5))
    S_adv = 2 * (maxU * dt)/(dx * (1-0.8))
    sysPrint("\t stability S_adv = %.06f, ideally less than 0.5" %S_adv)
    
    sysPrint('\t\tSeconds to run coupled step: %.4f' % (time.time() - start_time))

    # Now start all over again
