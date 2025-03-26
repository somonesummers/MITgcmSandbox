#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 12:39:00 2023

@author: jason
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 19:15:01 2023

@author: jason
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 14:04:35 2023

@author: jason
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
if(resetStart):
    input("Resetting this experiment directory, DELETING FILES AND FIGS. Confirm before continuing...")
    #Reset iceberg files
    # os.system('rm input/icebergs_length.bin input/icebergs_depths.bin input/icebergs_widths.bin')
    # os.system('cp input/icebergs_length_init.bin input/icebergs_length.bin')
    # os.system('cp input/icebergs_depths_init.bin input/icebergs_depths.bin')
    # os.system('cp input/icebergs_widths_init.bin input/icebergs_widths.bin')

    #Clean out itermediate save states and figures
    os.system("find couplingResults/MITgcmRun_*.pickle ! -name 'MITgcmRun_00000.pickle' -type f -exec rm {} +")
    os.system("rm couplingResults/*.txt")
    os.system("rm input/pickup*.data input/pickup*.meta")
    os.system("rm figs/*")

    #Reset data file
    for line in fileinput.input('input/data'):
            if "nIter0=" in line:
                startIter = float(line[8:-2])
            elif "endTime=" in line:
                endTime = float(line[9:-2])
    sysPrint('\tadjust start iteration %i to %i' %(int(startIter),int(0)))
    replaceAll('input/data','nIter0=%i' %(int(startIter)), 'nIter0=%i' % int(0))
    sysPrint('\tadjust end time %i to %i' %(int(endTime),int(86400)))
    replaceAll('input/data','endTime=%i' %(int(endTime)), 'endTime=%i' %(int(86400)))

    #Set melange to length consistent with initial glaciome1d size
    os.system('python advectBergs.py >> couplingResults/out.txt')

    #Run intial MITgcm, this resets the results folder
    os.system('bash ../makeRun.sh')
elif(freshStart):
    sysPrint('Running from fresh start...')
    input("Resetting the MITgcm directory. Confirm before continuing...")
    #Make couplingResults directory if not there already
    os.system("mkdir -p couplingResults")
    os.system('bash ../makeRun.sh')
else:
    sysPrint('Running from existing states...')


for ii in range(iterationsToRun):

    #Find time steps to take
    start_time = time.time()

    maxStep = 0
    for file in os.listdir('results'):
        # print(file)
        if "BRGFlx.0" in file:
            words = file.split(".")
            # print(words[1])  
            if int(words[1]) > maxStep:
                maxStep = int(words[1])

    dt = 0.0   
    for line in fileinput.input('input/data'):
            if "deltaT=" in line:
                dt = float(line[8:-2])
    # print('\tdt is loaded as', dt)


    x = np.squeeze(mds.rdmds("results/XC")[0,:])
    y = np.squeeze(mds.rdmds("results/YC")[:,0])
    dx = 2*x[0]
    dy = 2*y[0]
    ny = len(y)

    dataMITgcm = mds.rdmds("results/BRGFlx", maxStep)

            
    # n_pts = 21 # number of grid points
    # L = 10e3 # ice melange length
    # Ut = 6e3 # glacier terminus velocity [m/a]; treated as a constant
    # Uc = 6e3 # glacier calving rate [m/a]; treated as a constant
    # Ht = 800 # terminus thickness
    # n = 101 # number of time steps
    # dt = 0.002# 1/(n_pts-1)/10 # time step [a]; needs to be quite small for this to work
    # H0_manual = 100 #posit H0 start value, set to None to use default (75 m)

    # # specifying fjord geometry
    # X_fjord = np.linspace(0e3,400e3,101)
    # Wt = 4000
    # W_fjord = Wt + 0/10000*X_fjord
    # B = -1.0*constant.daysYear

    # first run to steady state
    # data = glaciome(n_pts, dt, L, Ut, Uc, Ht, B, X_fjord, W_fjord, H0 = H0_manual)
    #
    # OR load existing file as start
    files = sorted(glob.glob('couplingResults/MITgcmRun_[0-9][0-9][0-9][0-9][0-9].pickle'))
    # print(files)
    with open(files[-1], 'rb') as file:
        data = pickle.load(file)
        file.close()
        index = int(str(file).split(".")[1].split("_")[1])
        
    sysPrint('=== glaciome1D Running, load MITgcm data from end of day %.2f, using melange index %05i === ' %(int(maxStep)*dt/(24*3600),index))
    if(int(maxStep)*dt/(24*3600) - index != 1):
        raise Exception('MITgcm and glaciome1D are out of sync, glaciome index should be 1 behind MITgcm day')
    # os.system("echo ""glaciome1D Running, load MITgcm data from end of day %.2f, using melange index %05i"" > out.txt" %(int(maxStep)*dt/(24*3600),index))
    # data is of list ['BRGfwFlx','BRGhtFlx','BRGmltRt','BRG_TauX','BRG_TauY']
    # b_mitgcm = -1*np.nanmean(np.nansum(dataMITgcm[2,:,:,1:],axis=0),axis=0) # sum vertically, average across fjord
    b_mitgcm = -1*np.nanmean(np.nansum(dataMITgcm[0,:,:,1:],axis=0),axis=0)/(dx*dy*icebergCoverLambda)*(24*3600)*1000/917
    x_mitgcm = x[1:]-x[1] #ensure starts at 0, MITgcm has a glacier for first cell

    meltHelper = b_mitgcm.copy()
    meltHelper[meltHelper==0] = np.nan
    # meltHelper2 = b_mitgcm2.copy()
    # meltHelper2[meltHelper2==0] = np.nan
    sysPrint('\tMelt Rates min/mean/max: %.2f, %.2f, %.2f [m/day]:' %(np.nanmin(meltHelper),np.nanmean(meltHelper),np.nanmax(meltHelper)))
    # sysPrint('\tMelt Rates2 min/mean/max: %.2f, %.2f, %.2f[m/day]:' %(np.nanmin(meltHelper2),np.nanmean(meltHelper2),np.nanmax(meltHelper2)))
    # plt.figure()
    # plt.plot(x_mitgcm,b_mitgcm)
    # plt.show()
    # plt.close()
    data.dt = 1.0/365 # 1 day
    data.X_externalGrid = x_mitgcm
    data.B_externalGrid = b_mitgcm*365 #days to years
    sysPrint('\tPrevious Length: %.3f, H0: %.3f'%(data.L,data.H0))
    data.prognostic(method='lm')
    sysPrint('\tNew      Length: %.3f, H0: %.3f'%(data.L,data.H0))
    sysPrint('\t\tSeconds to run GLACIOME step: %.4f' % (time.time() - start_time))
    H = np.concatenate(([data.H0], data.H, [data.HL]))
    X = data.X
    X_ = np.concatenate(([data.X[0]],data.X_,[data.X[-1]]))
    np.save('input/melangeH',H)
    np.save('input/melangeX',X_)

    data.save('couplingResults/MITgcmRun_%05i.pickle' %(index+1))
    sysPrint('\tSaved as MITgcmRun_%05i.pickle' %(index+1))

    sysPrint('\tPrepping MITgcm for day %.2f' %(index+2))
    newStartTime = (index+1)*24*3600
    oldStartTime = (index)*24*3600
    sysPrint('\tmoving pickup')
    shutil.copy2('results/pickup.%010i.001.001.data' %int(newStartTime/dt), 'input')
    shutil.copy2('results/pickup.%010i.001.001.meta' %int(newStartTime/dt), 'input')

    sysPrint('\tadjust start iteration %i to %i' %(int(oldStartTime/dt),int(newStartTime/dt)))
    replaceAll('input/data','nIter0=%i' %(int(oldStartTime/dt)), 'nIter0=%i' % int(newStartTime/dt))
    sysPrint('\tadjust end time %i to %i' %(int(newStartTime),int(newStartTime+24*3600)))
    replaceAll('input/data','endTime=%i' %(int(newStartTime)), 'endTime=%i' %(int(newStartTime+24*3600)))

    os.system('python advectBergs.py >> couplingResults/out.txt')

    os.chdir("results")
    os.system('./mitgcmuv >> ../couplingResults/OutMITgcm%05i.txt' %(index+1))
    os.chdir("../")
    sysPrint('\t\tSeconds to run coupled step: %.4f' % (time.time() - start_time))
