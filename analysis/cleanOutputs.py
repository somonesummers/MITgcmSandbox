from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import xarray as xr

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 300
cleanPNGs = True

folders = ['.']
colors = ['xkcd:blue','xkcd:green']
resultFolder = '/results'
fileEnding = ""

if(os.path.isfile('input/plotHelperLocal.py')):
    sys.path.append('input')
    from plotHelperLocal import *
    print('Found experiment plotting settings')
elif(os.path.isfile('../plotHelper.py')):
    sys.path.append('../')
    print('no custom plotting settings, using local default')
    from plotHelper import *
else:  
    print('no defaults found')
print('Plot DPI:',plotDPI,'; clean PNGs?',cleanPNGs)



#Total melt over time
for j in range(len(folders)):
    folder = folders[j]

    dt = 0.0   
    for line in fileinput.input('%s/input/data' %folder):
            if "deltaT=" in line:
                dt = float(line[8:-2])
    print('dt is loaded as', dt)

    dz = np.load('input/dz.npy')
    x = mds.rdmds("results/XC")
    y = mds.rdmds("results/YC")
    z =mds.rdmds("results/RC")
    hFacC = mds.rdmds("results/hFacC")

    mds.wrmds('results/XC',x,itr=None, dataprec='float32')
    mds.wrmds('results/YC',y,itr=None, dataprec='float32')
    mds.wrmds('results/RC',z,itr=None, dataprec='float32')
    mds.wrmds('results/XC',x,itr=None, dataprec='float32')

    z = np.squeeze(z)
    #Find time steps to take
    maxStep = 0
    sizeStep = 1e10
    startStep = 1e10

    for file in os.listdir('%s%s' %(folder,resultFolder)):
        # print(file)
        if "BRGFlx.0" in file:
            words = file.split(".")
            # print(words[1])  
            if int(words[1]) > maxStep:
                maxStep = int(words[1])
            if int(words[1]) < startStep and int(words[1]) > 0:
                startStep = int(words[1])
            if abs(int(words[1]) - startStep) < sizeStep and abs(int(words[1]) - startStep) > 0:
                sizeStep = abs(int(words[1]) - startStep)
    
    maxFrames = 50
    if((maxStep-startStep)/sizeStep > maxFrames):   #if more than # frames, downscale to be less than #
        dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/maxFrames)
        print('Reducing time resolution by', dwnScale)
        sizeStep = sizeStep * dwnScale
    
    print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)




    dynName = ['BRGFlx','dynDiag','ptraceDiag']

    timeSteps = np.arange(startStep, maxStep + 1, sizeStep)
    for k in range(len(dynName)):
        print(dynName[k])
        for i in range(len(timeSteps)):
            data = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[k]), timeSteps[i])
            mds.wrmds('results/%s' %dynName[k],data,itr=timeSteps[i], dataprec='float32')
        

