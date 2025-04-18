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

MANUAL_PLUME_RATE = 25 #m^3/s

for j in range(len(folders)):
    folder = folders[j]

    dt = 0.0   
    for line in fileinput.input('%s/input/data' %folder):
            if "deltaT=" in line:
                dt = float(line[8:-2])
    print('dt is loaded as', dt)
    dz = np.load('input/dz.npy')

    hFacC = mds.rdmds("results/hFacC")

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
    dynName = ['BRGFlx']
    name = ['BRGfwFlx']
    units = ["[m^3/s]"]

    timeSteps = np.arange(startStep, maxStep + 1, sizeStep)
    fwOverTime = np.zeros(np.shape(timeSteps))
    plumeOverTime = np.zeros(np.shape(timeSteps))
    TracerPlumeTime = np.zeros(np.shape(timeSteps))
    TracerBergTime = np.zeros(np.shape(timeSteps))
    for i in range(len(timeSteps)):
        data = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), timeSteps[i])
        if i == 0:
        	fwOverTime[i] = np.nansum(data[0,:,:,:]) * sizeStep * dt
        else:
	        fwOverTime[i] = fwOverTime[i-1] + np.nansum(data[0,:,:,:]) * sizeStep * dt
        plumeOverTime[i] = MANUAL_PLUME_RATE*(i)*sizeStep*dt
        dataTracer = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'ptraceDiag'), timeSteps[i])
        dataTracer = dataTracer * dz[None,:,None,None] * 400 * 400 * hFacC
        TracerPlumeTime[i] = np.nansum(dataTracer[0,:,:,:])
        TracerBergTime[i] = np.nansum(dataTracer[1,:,:,:])

    plt.plot(timeSteps*dt/86400,fwOverTime,label='Berg Melt',color='xkcd:blue')
    plt.plot(timeSteps*dt/86400,plumeOverTime,label='Plume Water (manually set)',color='xkcd:red')
    plt.plot(timeSteps*dt/86400,TracerBergTime,label='Berg Trace',color='xkcd:blue',linestyle='--')
    plt.plot(timeSteps*dt/86400,TracerPlumeTime,label='Plume Trace',color='xkcd:red',linestyle='--')
plt.grid(alpha=.5)
plt.legend()
plt.title('Fresh Water Over Time')
plt.ylabel('Total water [m^3]')
plt.xlabel('Time [days]')
plt.savefig('figs/waterMass.png', format='png',dpi=plotDPI)
plt.show()
plt.close()
