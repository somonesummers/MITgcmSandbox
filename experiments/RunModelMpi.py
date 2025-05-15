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

# Handler for time limit on glaciome1d step
def handler(signum, frame):
    sysPrint("WARNING glaciome1d seems frozen, sad. Revert to old state")
    raise Exception("Timeout")

# Track exceptions
glmeWarningFlag = False
glmeWarningCount = 0
mitgcmWarningFlag = False
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
    sysPrint('\tadjust end time %i to %i' %(int(endTime),int(86400)))
    replaceAll('input/data','endTime=%i' %(int(endTime)), 'endTime=%i' %(int(86400)))

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
    endIter = int((24*3600)/dt)
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

    # Load mélange geometry

    files = sorted(glob.glob('couplingResults/MITgcmRun_[0-9][0-9][0-9][0-9][0-9].pickle'))
    # print(files)
    with open(files[-1], 'rb') as file:
        data = pickle.load(file)
        file.close()
        index = int(str(file).split(".")[1].split("_")[1])
        
    # make sure nothing funny is going on. Depending on when last cycle failed you may have to delete the highest numbered .pickle file
    if(int(maxStep)*dt/(24*3600) - index != 1):
        if(int(maxStep)*dt/(24*3600) - index == 0):
            #sometimes on restart MITgcm is forward, so autodelete and retake glaciome1d step if so
            sysPrint('glaciome1d 1 step ahead, deleting couplingResults/MITgcmRun_%05i.pickle' %index)
            os.system('rm couplingResults/MITgcmRun_%05i.pickle' %index)
            files = sorted(glob.glob('couplingResults/MITgcmRun_[0-9][0-9][0-9][0-9][0-9].pickle'))
            # print(files)
            with open(files[-1], 'rb') as file:
                data = pickle.load(file)
                file.close()
                index = int(str(file).split(".")[1].split("_")[1])
        else:
            sysPrint('ERROR load MITgcm data from end of day %.2f, using melange index %05i' %(int(maxStep)*dt/(24*3600),index))
            raise Exception('MITgcm (BRGFlx_) and GLACIOME (MITgcmRun_)are out of sync, GLACIOME index should be 1 behind MITgcm day')
    
    dataCopy = copy.deepcopy(data) #grab copy just in case we revert
    sysPrint('=== GLACIOME Running, load MITgcm data from end of day %.2f, using melange index %05i === ' %(int(maxStep)*dt/(24*3600),index))
    # Calculate melt rate. This is the net Freshwaterflux (sum verically, avg across)/ surface area of melange
    # It is very imporant to extend meltrate beyond mélange with value of last cell. Otherwise mélange will grow
    # continuiously. 

    lambdaHelper = 1-openFrac[0,1:-1,1:]
    lambdaHelper[lambdaHelper == 0] = 0.01
    # data is of list ['BRGfwFlx','BRGhtFlx','BRGmltRt','BRG_TauX','BRG_TauY']
    # Vertically sum all of the fresh water fluxes [m^3/s] cells/ divide by surface area 
    # (lamdba helper and dx*dy) of bergs in that cell [m/s]. dx*dy is outside the sum 
    # because it is constant, but still in the demnominator. Now a 'effective vertical melt rate' over surface
    # area at surface. Then average this value across the fjord, and multiply by seconds in a day to get 
    # vertical melt per day, the value glaciome wants. 
    b_mitgcm = -1*np.nanmean(np.nansum(dataMITgcm[0,:,1:-1,1:],axis=0)/lambdaHelper,axis=0)/(dx*dy)*(24*3600)

    if(np.mean(b_mitgcm)==0):
        # This means MITgcm has had an overflow. We can try once to re-advect bergs and try again. To do this, we delete the 
        # last BRGFlx file and skipping the rest of this loop with a CONTINUE statement. This consumes one iteration of the 
        # main loop that we will not get back. Also increment counter of how many of these weve had to ensure we dont get more
        # than 1 in a row
        if(mitgcmWarningFlag): #throw error if already in erroring state
            sysPrint('ERROR MITgcm failed 2x in a row')
            raise Exception('ERROR MITgcm failed 2x in a row')
        mitgcmWarningFlag = True #This is our first error in a row, activate error flag, count and warn
        mitgcmWarningCount = mitgcmWarningCount + 1
        sysPrint('\t\t WARNING MITgcm NANs, shuffling bergs and trying again. Failure %i' %mitgcmWarningCount)
        os.system('rm results/BRGFlx.%010i.*' %maxStep) #Delete last BRGFlx
        continue # This skips rest of loop. Next loop will delete last GLACIOME automatically, 
                 # run GLACIOME again, advect bergs, run MITgcm and hopefully resolve error. This technically does
                 # iterate the Melange for 2 days in 1 day of Ocean. Use carefully. Shorter timesteps should reduce
                 # this error from occuring.
    else:
        mitgcmWarningFlag = False

    # Now error are dealt with, We fill in the tail of 0s with the value of the melange toe 
    b_mitgcm[b_mitgcm == 0] = b_mitgcm[b_mitgcm != 0][-1]# if MITgcm has overflow error, this line fails

    x_mitgcm = x[1:]-x[1] #ensure starts at 0, MITgcm has a glacier for first cell
    
    sysPrint('\tMelt Rates min/mean/max: %.2f, %.2f, %.2f [m/day]:' %(np.min(b_mitgcm),np.mean(b_mitgcm),np.max(b_mitgcm)))
    # plt.figure()
    # plt.plot(x_mitgcm,b_mitgcm)
    # plt.show()
    # plt.close()
    data.dt = 1.0/365 # 1 day step for glaciome1d, as we will take 1 step
    data.X_externalGrid = x_mitgcm
    data.B_externalGrid = b_mitgcm*365 #days to years

    # Ready to now run GLACIOME1D. We run 1 timestep with dt set to the coupling timestep. If this creates CLF issues
    # we can run more steps with a reduced dt, or couple with MITgcm more frequently, but this isn't likely to be needed. 
    oldL = data.L
    oldH = data.H0

    if(forceMelange):
        data.Uc = 6000 - 2 * index
        data.Ut = 6000 - 2 * index # Must change together unless terminus moving
        sysPrint('\tNew Uc: %.3f m/year'%(data.Uc))
    sysPrint('\tPrevious Length: %.3f, H0: %.3f'%(oldL, oldH))
    # data.steadystate()

    signal.signal(signal.SIGALRM,handler) #setting a timer for glaciome1d
    signal.alarm(600) # if its not done in 10 minutes, its probably at minimum size. 
                    # if this is happening for non-trivially small melange, something is off.
    try:
        data.prognostic(method='lm')
    except (Exception):
        data = dataCopy
        glmeWarningFlag = True
    signal.alarm(0) #turn the alarm off, because we made it
    sysPrint('\tNew      Length: %.3f, H0: %.3f, ∆L: %.3f, ∆H0: %.3f'%(data.L,data.H0,data.L-oldL,data.H0-oldH,))
    sysPrint('\t\tSeconds to run GLACIOME step: %.4f' % (time.time() - start_time))
    
    if(data.H0 < 25.0 or data.L < 300): #limit small size of melang, glaciome1d gets slow 
        # This is really a secondary limit for the same reason as the time limit above
        # likely, we should only have the time limit and not this one.
        data = dataCopy
        glmeWarningFlag = True
        sysPrint('\t\t WARNING melange below minimum size, revert to old size')
    
    # If we had a failure, increment counter. If succeeded, reset counter
    if(glmeWarningFlag):
        glmeWarningCount = glmeWarningCount + 1
        sysPrint('\t\t\t Successive GLACIOME failures: %i ' %glmeWarningCount)
        glmeWarningFlag = False
    else:
        glmeWarningCount = 0

    if(glmeWarningCount > 10 and data.L > 1000):
        sysPrint('ERROR GLACIOME failed to run 10 time in a row for nontrivial melange size')
        raise Exception('ERROR GLACIOME failed to run 10 time in a row for nontrivial melange size')


    # save itermediate steps for plotting
    data.save('couplingResults/MITgcmRun_%05i.pickle' %(index+1))
    sysPrint('\tSaved as MITgcmRun_%05i.pickle' %(index+1))

    # Preparing MITgcm for next run. 
    # We adust the data file, start and end fields. 

    sysPrint('\tPrepping MITgcm for day %.2f' %(index+2))
    newStartTime = (index+1)*24*3600
    oldStartTime = (oldStartIter)*dt #read directly from input/data now
    # we delete the old pickups to reduce file size and file count. We could strategically leave them every
    # 100 iterations or so, but currently we do not. 
    os.system('rm results/pick*.%010i.*' %int((oldStartTime-24*3600)/dt)) #save last one, but delete 2 ago

    sysPrint('\tadjust start iteration %i to %i' %(oldStartIter,int(newStartTime/dt)))
    replaceAll('input/data','nIter0=%i' %(int(oldStartIter)), 'nIter0=%i' % int(newStartTime/dt))
    sysPrint('\tadjust end time %i to %i' %(int(newStartTime),int(newStartTime+24*3600)))
    replaceAll('input/data','endTime=%i' %(int(newStartTime)), 'endTime=%i' %(int(newStartTime+24*3600)))

    if(index > 1 and index % 100 == 0):
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

    endIter = int((newStartTime + 24*3600)/dt)
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
