#!/usr/bin/env python3
"""
Created on Wed Mar 26 2025
Updated on Jun 2026 to be work on Tufts Pax and GATech PACE
@author: psummers8

Main Script for running mélange/MITgcm coupled model.
Depends on advectBergs.py, this could (should) be incorperated. 
Depends on glaciome1d class. This should stand alone.



"""

import numpy as np
from matplotlib import pyplot as plt
import sys
import fileinput
import shutil
from MITgcmutils import mds
import os
import re
# sys.path.append('/hdd/glaciome/models/glaciome1D')
# sys.path.insert(0, '')
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
sys.path.append('/cluster/home/psumme03/glaciome1d')
from glaciome1D import constants, glaciome
import glob
import pickle
import time
import copy
import signal
import platform


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

OSX = platform.system()
current_directory = os.getcwd()
baseDir = ''
CondaDir = ''
if OSX == 'Darwin':
    CondaDir = '~/anaconda3'
elif "cluster" in current_directory:
    CondaDir = '~/.conda'
    OSX = 'Tufts'
elif "hcoda" in current_directory:
    CondaDir = '~/.conda'
    OSX = 'PACE'
else:
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

# Track exceptions
glmeWarningFlag = False
glmeWarningCount = 0
mitgcmWarningFlag = 1 #if it fails first run, should error
mitgcmWarningCount = 0


if(resetStart):
    raise Exception("Restart start no longer supported")
elif(freshStart): #This distinguises between a new coupled run, or continuing a paused coupled run
    #Make couplingResults directory if not there already
    os.system("mkdir -p couplingResults")
    os.system("rm -f couplingResults/out.txt")
    sysPrint('Running from fresh start...')
    sysPrint("Assuming fresh MITgcm directory. Steps to run: %i..." %iterationsToRun)
    time.sleep(1) # this is a brief pause that lets you panic kill the job as the next steps nukes the results folder

    #mpirun flags, these help reduce minor errors and unneeded output
    mpi_flags = ''
    if(OSX == 'Tufts'):
        mpi_flags = '--quiet --mca btl_vader_single_copy_mechanism none'
    elif(OSX == 'PACE'):
        mpi_flags = '-v'
    elif(OSX == 'Darwin'):
        mpi_flags = '-v'
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

    prefixes = ['BRGFlx','dynDiag','ptraceDiag','plumeDiag']#,'presDiag','momDiag','heatDiag']
    # prefixes = ['presDiag','momDiag','heatDiag']
    endIter = int((24*3600)/dt)
    # sysPrint('\tcondensing diagnostic tile files to global files iter:%i' %endIter)
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

    # Load mélange geometry

    files = sorted(glob.glob('couplingResults/MITgcmRun_[0-9][0-9][0-9][0-9][0-9].pickle'))
    # print(files)
    if(len(files) == 0):
        raise Exception('No mélange files found')
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
    lambdaHelper[lambdaHelper == 0] = 0.01 #set minimum for icefrac to smooth edge cases, maybe un-needed?
    # data is of list ['BRGfwFlx','BRGhtFlx','BRGmltRt','BRG_TauX','BRG_TauY']
    # Vertically sum all of the fresh water fluxes [m^3/s] cells/ divide by surface area 
    # (lamdba helper and dx*dy) of bergs in that cell [m/s]. dx*dy is outside the sum 
    # because it is constant, but still in the demnominator. Now a 'effective vertical melt rate' over surface
    # area at surface. Then average this value across the fjord, and multiply by seconds in a day to get 
    # vertical melt per day, the value glaciome wants. 
    b_mitgcm = -1*np.nanmean(np.nansum(dataMITgcm[0,:,1:-1,1:],axis=0)/lambdaHelper,axis=0)/(dx*dy)*(24*3600)

    if(np.mean(b_mitgcm) > -1e-4 or np.mean(b_mitgcm) < -1000): #sometimes returns 0 or -inf both are problems
        sysPrint(f'\t\tMITgcm melt rates looks problematic min/ave/min: {np.min(b_mitgcm):.2f}/{np.mean(b_mitgcm):.2f}/{np.max(b_mitgcm):.2f}')
        # This means MITgcm has had an overflow. We can try once to re-advect bergs and try again. To do this, we delete the 
        # last BRGFlx file and skipping the rest of this loop with a CONTINUE statement. This consumes one iteration of the 
        # main loop that we will not get back. Also increment counter of how many of these weve had to ensure we dont get more
        # than 1 in a row
        mitgcmWarningFlag = mitgcmWarningFlag + 2 #increment by 2, this is the error window
        if(mitgcmWarningFlag > 2): #throw error if already in erroring state
            sysPrint('ERROR MITgcm failed 2x in a row')
            raise Exception('ERROR MITgcm failed 2x in a row')
        mitgcmWarningCount = mitgcmWarningCount + 1
        sysPrint('\t\t WARNING MITgcm NANs, shuffling bergs and trying again. Failure %i' %mitgcmWarningCount)
        os.system('rm results/BRGFlx.%010i.*' %maxStep) #Delete last BRGFlx
        #reset ICEBERG mitgcm files to previous state
        for tmp_name in ['bergMask','numBergsPerCell','openFrac','totalBergArea','meltMask','driftMask','barrierMask','icebergs_depths','icebergs_widths','icebergs_length']:
            os.system(f'cp -v input/temp_{tmp_name}.bin input/{tmp_name}.bin' )
        continue # This skips rest of loop. Next loop will delete last GLACIOME automatically, 
                 # run GLACIOME again, advect bergs, run MITgcm and hopefully resolve error. This technically does
                 # iterate the Melange for 2 days in 1 day of Ocean. Use carefully. Shorter MITgcm timesteps (dt) should reduce
                 # this error from occuring.
    else:
        mitgcmWarningFlag = 0 if mitgcmWarningFlag == 0 else mitgcmWarningFlag - 1

    # Now error are dealt with, We fill in the tail of 0s with the value of the melange toe 
    if(False): #option to bias mélange to grow by having no melt beyond toe of mélange
        sysPrint('\t == forcing growth == ')
        b_mitgcm[b_mitgcm != 0] = 0
    else:
        b_mitgcm[b_mitgcm == 0] = b_mitgcm[b_mitgcm != 0][-1]# if MITgcm has overflow error, this line fails
    x_mitgcm = x[1:]-x[1] #ensure starts at 0, and ignore first cell. MITgcm has a glacier for first cell
    
    sysPrint('\tMelt Rates min/mean/max: %.2f, %.2f, %.2f [m/day]:' %(np.min(b_mitgcm),np.mean(b_mitgcm),np.max(b_mitgcm)))
    # plt.figure()
    # plt.plot(x_mitgcm,b_mitgcm)
    # plt.show()
    # plt.close()
    data.dt = 1.0/365 # 1 day step for glaciome1d, as we will take 1 step
    data.X_externalGrid = x_mitgcm + data.X[0] # if terminus moves in GLACIOME, the melt grid must too
    data.B_externalGrid = b_mitgcm*365 #days to years

    # Ready to now run GLACIOME1D. We run 1 timestep with dt set to the coupling timestep. If this creates CLF issues
    # we can run more steps with a reduced dt, or couple with MITgcm more frequently, but this isn't likely to be needed. 
    oldL = data.L
    oldH = data.H0

    if(forceMelange):
        alpha = couplingAlpha
        beta = couplingBeta
        U0 = couplingUc0
        F = data.force()
        ucRamp = 0
       # if(index > 35):
        #    ucRamp = (index - 35) * 5
         #   sysPrint(f'\tRamping on: {ucRamp}')
        data.Uc = U0 + ucRamp - alpha * F - beta * data.X[0]
        data.Ht = 600 - data.X[0] * beta
    sysPrint('\tPrevious Length: %.3f, H0: %.3f'%(oldL, oldH))
    # data.steadystate()

    signal.signal(signal.SIGALRM,handler) #setting a timer for glaciome1d
    signal.alarm(6000) # if its not done in 100 minutes, its probably at minimum size. 
                    # if this is happening for non-trivially small melange, something is off.
                    # Can happen when mélange thickness increases down fjord, caused by freezing which should be fixed now 
    try:
        data.prognostic(method='lm')
    except Exception as e:
        sysPrint(str(e))
        data = dataCopy
        data.t = data.t + data.dt #step forward one day to keep up with MITgcm
        glmeWarningFlag = True
    signal.alarm(0) #turn the alarm off, because we made it
    sysPrint('\tNew      Length: %.3f, H0: %.3f, ∆L: %.3f, ∆H0: %.3f'%(data.L,data.H0,data.L-oldL,data.H0-oldH,))
    sysPrint(f"\tterm loc/depth {data.X[0]:5.0f}/{data.Ht:5.1f} m, calving/glacier speed: {data.Uc:0.0f}/{data.Ut:0.0f} m/yr")
    sysPrint('\t\tSeconds to run GLACIOME step: %.4f' % (time.time() - start_time))
    
    if(data.L < 300): #limit small size of melang, glaciome1d gets slow 
        # This is really a secondary limit for the same reason as the time limit above
        # likely, we should only have the time limit and not this one.
        data = dataCopy
        data.t = data.t + data.dt #step forward one day to keep up with MITgcm
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
    os.system('rm results/pick*.%010i.*' %int((oldStartTime-48*3600)/dt)) #save last two, but delete 3 ago
    # we save 2 as we let MITgcm fail 2x in a row before failing the whole program, so this lets us have
    # an old good state (hopefully)

    sysPrint('\tadjust start iteration %i to %i' %(oldStartIter,int(newStartTime/dt)))
    replaceAll('input/data','nIter0=%i' %(int(oldStartIter)), 'nIter0=%i' % int(newStartTime/dt))
    sysPrint('\tadjust end time %i to %i' %(int(newStartTime),int(newStartTime+24*3600)))
    replaceAll('input/data','endTime=%i' %(int(newStartTime)), 'endTime=%i' %(int(newStartTime+24*3600)))

    if(False and index > 1 and index % 100 == 0): 
    #the FALSE here is because this not done anymore, but you can flip it if you would like
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
    os.system(f"{CondaDir}/envs/MITgcm/bin/python advectBergsCombo.py >> couplingResults/out.txt")

    os.chdir("results")
    maxTry = 3
    # Sometime PACE specifically will timeout when queueing SRUN jobs, this kills the whole run and is a bummer,
    # so here we allow the script to try a few times to get this to go. 
    for attempt_i in range(maxTry):
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
                sys_ret = os.system(f"mpirun {mpi_flags} -n {nPy*nPx} ./mitgcmuv")
            elif(OSX == 'PACE'): 
                sys_ret = os.system(f"srun {mpi_flags} ./mitgcmuv") #srun is smart and allocated cores automatically
            elif(OSX == 'Darwin'):
                sys_ret = os.system(f"mpirun {mpi_flags} -n {nPy*nPx} ./mitgcmuv")
            else:
                raise ExceptionType("unknown OSX or running location, please configure")
        else:
            sys_ret = os.system(f"./mitgcmuv > ../MITgcm_output.txt")
        sysPrint(f"MITgcm run return is {sys_ret}")
        if(sys_ret == 0): #MITgcm ran successfully, so lets continue. Specifically sys_ret = 256 seems to be the 
                          # srun: error: Unable to confirm allocation for job 9983843: Socket timed out on send/recv operation
                          # srun: Check SLURM_JOB_ID environment variable. Expired or invalid job 9983843
            break
    os.chdir("../")
    
    # Create global files to reduce file counts
    # MITgcm MPI saves outputs per tile, per timestep. For 20 cores this makes 40 files per diagnostic
    # per timestep, which quickly become insane for long runs. This step collects all those tile files into
    # 2 global files (*.data, *.meta), then deletes the tiles files. This reduces files count by a factor of 20
    # This VASTLY improves data transfer and compression speeds. 
    prefixes = ['BRGFlx','dynDiag','ptraceDiag','plumeDiag']#,'presDiag','momDiag','heatDiag']

    endIter = int((newStartTime + 24*3600)/dt)
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
