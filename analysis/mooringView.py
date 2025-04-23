from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import xarray as xr
import gsw
import argparse

parser = argparse.ArgumentParser(description='Plot mooring view at xCrossSection')
parser.add_argument('xCrossSection', nargs='?', const=0.0, type=float,
                    help='optional x location [m]')
args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 300
cleanPNGs = True

folder = '.'
resultFolder = '/results'

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

if(args.xCrossSection != None):
    print('** Manual xCrossSection detected **')
    xCrossSection = args.xCrossSection

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
nz = len(z)
hFacC = mds.rdmds("results/hFacC")

xSlice = np.argmin(np.abs(x[0,:] - xCrossSection))
print('cross section is x =', x[0,xSlice],'index', xSlice)

#Find time steps to take
maxStep = 0
sizeStep = 1e10
startStep = 1e10

for file in os.listdir('%s%s' %(folder,resultFolder)):
    # print(file)
    if "dynDiag.0" in file:
        words = file.split(".")
        # print(words[1])  
        if int(words[1]) > maxStep:
            maxStep = int(words[1])
        if int(words[1]) < startStep and int(words[1]) > 0:
            startStep = int(words[1])
        if abs(int(words[1]) - startStep) < sizeStep and abs(int(words[1]) - startStep) > 0:
            sizeStep = abs(int(words[1]) - startStep)

if((maxStep-startStep)/sizeStep > 60):   #if more than 50 frames, downscale to be less than 50
    dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/60)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

timeSteps = np.arange(startStep, maxStep + 1, sizeStep)
dynTime = np.zeros([5,nz,len(timeSteps)])

for i in range(len(timeSteps)):
    # ["Temp", "Sal", "U", "W", "V"]
    dataOcean = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'dynDiag'), timeSteps[i])
    dynTime[:,:,i] = np.nanmean(dataOcean[:,:,1:-1,xSlice],axis=2)

name = ["Temp", "Sal", "U"]
cbarLabel = ["[C]", "[PSU]", "[m/s]"]  
lvlList = [tempRange,saltRange,uRange]
cmList = [tempCmap,saltCmap,uCmap]

fig, axes = plt.subplots(3, 1, figsize=(12, 8), layout="constrained")
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]


for k in range(len(name)):
    lvl = lvlList[k]
    cm = cmList[k]
    cp = axes[k].contourf(timeSteps*dt/86400,
                    z,
                    dynTime[k,:,:],
                    lvl,
                    extend="both",
                    cmap=cm)
    cbar = plt.colorbar(cp)
    cbar.set_label(name[k] + ' ' + cbarLabel[k])
    if(k == len(name)-1):
        axes[k].set_xlabel('Time [days]')
    axes[k].set_xticks(range(0, int(timeSteps[-1]*dt/86400), 91))
    axes[k].set_ylabel('Depth [m]')
    axes[k].grid(alpha=.2)
plt.suptitle('Mooring at %.02f meters along fjord' %x[0,xSlice])
plt.savefig('figs/mooring%05i.png' %xCrossSection, format='png',dpi=plotDPI)
plt.show()
plt.close()

