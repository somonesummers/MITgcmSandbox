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
    dy = y[0,0]/2
    dx = x[0,0]/2
    z = np.squeeze(mds.rdmds("results/RC"))
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
    tempOverTime = np.zeros(np.shape(timeSteps))
    tempOverTime2 = np.zeros(np.shape(timeSteps))
    for i in range(len(timeSteps)):
        #name = ['BRGfwFlx','BRGhtFlx','BRGmltRt','BRG_TauX','BRG_TauY']
        data = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), timeSteps[i])
        #name = ["W", "Temp", "Sal", "CellMeltRate", "RadiusThickness","Flux","Density"]
        dataPlume = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'plumeDiag'), timeSteps[i])
        # ["Temp", "Sal", "U", "W", "V"]
        dataOcean = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'dynDiag'), timeSteps[i])
        dataOcean[0,data[0,:,:,:] == 0] = np.nan #average over just the melange

        dataOceanBelow = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'dynDiag'), timeSteps[i])
        dataOceanBelow[:,z > -300,:,:] = np.nan #nan above 300 meters depth
        dataOceanBelow[:,:,:,int(10000/dx):] = np.nan #beyond 10km in X to nan

        tempOverTime[i] = np.nanmean(dataOcean[0,:,:,:])
        tempOverTime2[i] = np.nanmean(dataOceanBelow[0,:,:,:])
        fwOverTime[i] = np.nansum(data[0,:,:,:])
        plumeOverTime[i] = np.nanmax(dataPlume[0,-1,:,:] * dataPlume[4,-1,:,:] * dy) #assume sheet here 
        # this is a little high as some water has entrained

    plt.plot(timeSteps*dt/86400,fwOverTime,label='Iceberg Meltwater',color='xkcd:blue')
    plt.plot(timeSteps*dt/86400,plumeOverTime,label='Plume Water',color='xkcd:red')
    plt.grid(alpha=.5)
    ax1 = plt.gca()
    ax2=ax1.twinx()
    ax2.grid(alpha=.1,color="xkcd:pumpkin")
    ax2.plot(timeSteps*dt/86400,tempOverTime,label='Mélange Avg Temp',color='xkcd:pumpkin')
    ax2.plot(timeSteps*dt/86400,tempOverTime2,label='Below 300 m Avg Temp',color='xkcd:twilight')
plt.grid(alpha=.5)
ax1.legend()
ax2.legend()
plt.title('Fresh Water Over Time')
ax1.set_ylabel('Total water [m^3]')
ax2.set_ylabel('Mélange Temperature [C]',color='xkcd:orange')
plt.xlabel('Time [days]')
plt.savefig('figs/waterFlux.png', format='png',dpi=plotDPI)
plt.show()
plt.close()

