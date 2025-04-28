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

    prefixes = ['Depth','DXC','DXF','DXG','DXV','DYC','DYF','DYG','DYU','hFacC','hFacS','hFacW',
                        'maskInC','maskInS','maskInW','RAC','RAS','RAW','RAZ','XC','XG','YC','YG']
    print('\tcondensing grid tile files to global files')
    for k in range(len(prefixes)):
        print('\t\t== %s ==' %prefixes[k])
        if(os.path.isfile('results/%s.001.001.meta'%prefixes[k])):
            dataTemp = mds.rdmds("results/%s"%(prefixes[k]))
            mds.wrmds('results/%s' %prefixes[k],dataTemp,dataprec='float32')
            os.system('rm results/%s.0*.0*' %(prefixes[k])) #picks out tile level files
        else:
            print('\t\t\t Already Clean')

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
    
    print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)


    dynName = ['BRGFlx','dynDiag','ptraceDiag','plumeDiag']
    print('\tcondensing diagnostic tile files to global files')
    timeSteps = np.arange(startStep, maxStep + 1, sizeStep)
    for k in range(len(dynName)):
        print(dynName[k])
        for i in range(len(timeSteps)):
            if(os.path.isfile('results/%s.%010i.001.001.meta'%(dynName[k],timeSteps[i]))):
                data = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[k]), timeSteps[i])
                mds.wrmds('results/%s' %dynName[k],data,itr=timeSteps[i], dataprec='float32')
                os.system('rm results/%s.%010i.0*.0*' %(dynName[k],timeSteps[i]))
        

