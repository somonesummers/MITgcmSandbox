from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import argparse

# Take input options. Some defaults are set here, so be aware
parser = argparse.ArgumentParser(description='Plot options for fresh water flux over time')
parser.add_argument('-s','--silent', action='count', default=0,
                    help='Option to silence showing of plots')
parser.add_argument('-n','--numFrames', nargs=1, default=[150],type=int,
                    help='Max number of samples from time series [defaut = 150]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('-sgd','--sgd', action='count', default=0,
                    help='optional specification to normalize by (SGD/1000), T,S excluded [default = off]')
parser.add_argument('-x','--xCrossSection', nargs=1, type=float,default = None,
                    help='optional slice location [m]')
parser.add_argument('-z','--zCrossSection', nargs=1, type=float,default = None,
                    help='optional slice location [m]')
args = parser.parse_args()


# Pick cross section to view from file or default
plotDPI = 150
manualMax = None

folders = ['.']
labels = ['glaciome']
colors = ['xkcd:blue','xkcd:gray','xkcd:red']
resultFolder = '/results'
fileEnding = ""

if(args.timeRange != None):
    fileEnding = f'{fileEnding}_T{args.timeRange[0]}_{args.timeRange[1]}' 
thresholds = [0.04,.04,.04]

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
print('Plot DPI:',plotDPI)



dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

y = np.squeeze(mds.rdmds("results/YC")[:,0])
x = np.squeeze(mds.rdmds("results/XC")[1])
z = np.squeeze(mds.rdmds("results/RC"))
dz = np.load("results/dz.npy")
ny = len(y)

if(args.zCrossSection == None):
    zDepth = np.min(z) #only look above this
else:
    zDepth = args.zCrossSection[0]
    fileEnding = fileEnding + f'_z{zDepth}'
zSlice = np.argmin(np.abs(z[:]- zDepth))
print('depth is z =', z[zSlice], 'index', zSlice)

if(args.xCrossSection == None):
    xReach = 0 #only look above this
else:
    xReach = args.xCrossSection[0]
    fileEnding = fileEnding + f'_x{xReach}'
xSlice = np.argmin(np.abs(x[:]- xReach))
print('reach is x =', x[xSlice], 'index', xSlice)


#extracting SGD
sgd = (args.sgd > 0)
if(sgd):
    print('Normalizing by SGD')
    fileEnding = fileEnding + 'SGDNorm'
    from scipy import interpolate
    from MITgcmutils import mds
    y = np.squeeze(mds.rdmds("results/YC")[:,0])
    dy = 2*y[0]
    for line in fileinput.input('input/data'):
        if "ExternForcingPeriod=" in line:
            period = float(line[21:-2])
        elif "ExternForcingCycle=" in line:
            cycle = float(line[20:-2])
    nt = int(cycle/period)
    plumeMask = np.fromfile('input/plumeMask.bin', dtype='>f8') #(2 is line, 3 is semi-cone)
    rad = np.fromfile('input/runoffRad.bin', dtype='>f8')
    rad = rad.reshape((nt,len(plumeMask)))
    vel = np.fromfile('input/runoffVel.bin', dtype='>f8')
    vel = vel.reshape((nt,len(plumeMask)))

    plumeType = np.max(plumeMask)
    cleanRad = rad[:,plumeMask == plumeType]
    cleanVel = vel[:,plumeMask == plumeType]
    seasonTime = np.linspace(0,cycle/86400,nt)
    if(plumeType == 2):
        runoff = dy*cleanRad[:,0]*cleanVel[:,0]
    else: #semi-rad
        runoff = np.pi*cleanRad**2*cleanVel
    f = interpolate.interp1d(seasonTime, runoff,fill_value='0')

figLabels = ["(a)","(b)","(c)","(d)","(d)","(f)","(g)","(h)"]

for j in range(len(folders)):
    folder = folders[j]

    #Find time steps to take
    maxStep = 0
    sizeStep = 1e10
    startStep = 1e10

    for file in os.listdir('%s%s' %(folder,resultFolder)):
        # print(file)
        if "momDiag.0" in file:
            words = file.split(".")
            # print(words[1])  
            if int(words[1]) > maxStep:
                maxStep = int(words[1])
            if int(words[1]) < startStep and int(words[1]) > 0:
                startStep = int(words[1])
            if abs(int(words[1]) - startStep) < sizeStep and abs(int(words[1]) - startStep) > 0:
                sizeStep = abs(int(words[1]) - startStep)

    if(args.timeRange != None):
        startStep = args.timeRange[0] * 86400 / dt
        if(args.timeRange[1] != 0):
            maxStep = args.timeRange[1] * 86400 / dt

    maxFrames = args.numFrames[0]
    if((maxStep-startStep)/sizeStep > maxFrames):   #if more than # frames, downscale to be less than #
        dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/maxFrames)
        print('Reducing time resolution by', dwnScale)
        sizeStep = sizeStep * dwnScale

    print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

    dynName = ['momDiag']
    # name = ['BRGfwFlx']
    # units = ["[m^3/s]"]


    weightsTmp = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), startStep)
    wghts = weightsTmp[0,:,:,:]
    for i in range(len(dz)):
        wghts[i,:,:] = dz[i]

    timeSteps = np.arange(startStep, maxStep + 1, sizeStep)
    # Depth view
    uDissDepth   = np.zeros([len(z),len(timeSteps)])
    uAdvecDepth  = np.zeros([len(z),len(timeSteps)])
    uCoriDepth   = np.zeros([len(z),len(timeSteps)])
    udPhiDepth   = np.zeros([len(z),len(timeSteps)])
    uBergDepth   = np.zeros([len(z),len(timeSteps)])
    spdDepth     = np.zeros([len(z),len(timeSteps)])
    # Length view
    uDissLength  = np.zeros([len(x),len(timeSteps)])
    uAdvecLength = np.zeros([len(x),len(timeSteps)])
    uCoriLength  = np.zeros([len(x),len(timeSteps)])
    udPhiLength  = np.zeros([len(x),len(timeSteps)])
    uBergLength  = np.zeros([len(x),len(timeSteps)])
    vBergLength  = np.zeros([len(x),len(timeSteps)])
    spdLength    = np.zeros([len(x),len(timeSteps)])
    for i in range(len(timeSteps)):
        data = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), timeSteps[i])
        dataBrg = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'BRGFlx'), timeSteps[i])
        dataOcean = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'dynDiag'), timeSteps[i])
        melangeMask = dataBrg[0,:,:,:].copy()
        melangeMask[melangeMask != 0] = 1
        melangeMask[data[0,:,:,:] < -1] = 1 #attemp to recapture freezing limit ice
        melangeMask[:,np.nansum(data[0,:,:,:],axis=0) == 0] = 0 #attemp to recapture freezing limit ice
        # melangeMask[:,:,1] = 0 #exclude glacier face and plume
        data[:, melangeMask == 0] = np.nan #nan all zero melt cells
        dataBrg[:, melangeMask == 0] = np.nan #nan all zero melt cells
        dataOcean[:, melangeMask == 0] = np.nan #nan all zero melt cells
        if(sgd):
            sgdFactor = f(i*dt)/1000
        else:
            sgdFactor = 1
        #depth values
        uDissDepth[:,i]   = np.nanmean(data[0,:,:,xSlice:],axis=(1,2))/sgdFactor
        uAdvecDepth[:,i]  = np.nanmean(data[1,:,:,xSlice:],axis=(1,2))/sgdFactor
        uCoriDepth[:,i]   = np.nanmean(data[2,:,:,xSlice:],axis=(1,2))/sgdFactor
        udPhiDepth[:,i]   = np.nanmean(data[3,:,:,xSlice:],axis=(1,2))/sgdFactor
        uBergDepth[:,i]   = np.nanmean(dataBrg[3,:,:,xSlice:] / (1020 * dz[:,None,None]),axis=(1,2)) /sgdFactor
        spdDepth[:,i]     = np.nanmean(np.sqrt(dataOcean[2,:,:,xSlice:]**2 + dataOcean[3,:,:,xSlice:]**2 +
                                dataOcean[4,:,:,xSlice:]**2),axis=(1,2))/sgdFactor
        #length values
        uDissLength[:,i]   = np.nanmean(data[0,:zSlice,:,:],axis=(0,1))/sgdFactor
        uAdvecLength[:,i]  = np.nanmean(data[1,:zSlice,:,:],axis=(0,1))/sgdFactor
        uCoriLength[:,i]   = np.nanmean(data[2,:zSlice,:,:],axis=(0,1))/sgdFactor
        udPhiLength[:,i]   = np.nanmean(data[3,:zSlice,:,:],axis=(0,1))/sgdFactor
        uBergLength[:,i]   = np.nanmean(dataBrg[3,:zSlice,:,:] / (1020 * dz[:zSlice,None,None]),axis=(0,1))/sgdFactor
        spdLength[:,i]     = np.nanmean(np.sqrt(dataOcean[2,:zSlice,:,:]**2 + dataOcean[3,:zSlice,:,:]**2 +
                                dataOcean[4,:zSlice,:,:]**2),axis=(0,1))/sgdFactor

## Next figure
fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]
# ax5 = axes[2,0]
# ax6 = axes[2,1]

momColorMap = 'cmo.balance'
wMomColorMap = 'cmo.curl'
uMomRange = np.linspace(-2e-4,2e-4,31)
# u2MomRange = np.linspace(-2e-6,2e-6,31)
wMomRange = np.linspace(-1e-5,1e-5,31)

zMin=np.min(z[~np.isnan(np.nanmean(uDissDepth,axis=1))])
cp = ax1.contourf(timeSteps*dt/86400,z,uDissDepth - uBergDepth,
                uMomRange,
                cmap=momColorMap,
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('[m/s^s]')
ax1.grid(alpha=.5)
ax1.set_title('uDiss')
ax1.set_ylabel('Depth [m]')
ax1.set_xlabel('Time [days]')
ax1.set_ylim([zMin,0])

cp = ax2.contourf(timeSteps*dt/86400,z,udPhiDepth,
                uMomRange,
                cmap=momColorMap,
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('[m/s^s]')
ax2.grid(alpha=.5)
ax2.set_title('udPhiDepth')
ax2.set_ylabel('Depth [m]')
ax2.set_xlabel('Time [days]')
ax2.set_ylim([zMin,0])

cp = ax3.contourf(timeSteps*dt/86400,z,uBergDepth,
                uMomRange,
                cmap=momColorMap,
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('[m/s^s]')
ax3.grid(alpha=.5)
ax3.set_title('uBergDepth')
ax3.set_ylabel('Depth [m]')
ax3.set_xlabel('Time [days]')
ax3.set_ylim([zMin,0])

cp = ax4.contourf(timeSteps*dt/86400,z,spdDepth,
                np.linspace(0,.1,11),
                cmap='cmo.speed',
                extend="both")
cp2 = ax4.contour(timeSteps*dt/86400,z,spdDepth,
                [0.04],
                colors='red',
                linewidths=1,
                linestyles='--',
                alpha=.5)
cbar = plt.colorbar(cp)
cbar.set_label('[m/s]')
ax4.grid(alpha=.5)
ax4.set_title('Speed')
ax4.set_ylabel('Depth [m]')
ax4.set_xlabel('Time [days]')
ax4.set_ylim([zMin,0])

plt.suptitle(f'View Beyond {xReach} m')
plt.savefig('figs/momDepthView%s.png' %fileEnding, format='png',dpi=plotDPI)
if(args.silent == 0):
    plt.show()
plt.close()

## Next figure
fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]
# ax5 = axes[2,0]
# ax6 = axes[2,1]

xMax=np.max(x[~np.isnan(np.nanmean(uDissLength,axis=1))])
cp = ax1.contourf(timeSteps*dt/86400,x,uDissLength - uBergLength,
                uMomRange,
                cmap=momColorMap,
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('[m/s^2]')
ax1.grid(alpha=.5)
ax1.set_title('uDissLength')
ax1.set_ylabel('Length [m]')
ax1.set_xlabel('Time [days]')
ax1.set_ylim([0,xMax])

cp = ax2.contourf(timeSteps*dt/86400,x,udPhiLength,
                uMomRange,
                cmap=momColorMap,
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('[m/s^2]')
ax2.grid(alpha=.5)
ax2.set_title('udPhiLength')
ax2.set_ylabel('Length [m]')
ax2.set_xlabel('Time [days]')
ax2.set_ylim([0,xMax])

cp = ax3.contourf(timeSteps*dt/86400,x,uBergLength,
                uMomRange,
                cmap=momColorMap,
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('[m/s^2]')
ax3.grid(alpha=.5)
ax3.set_title('uBergDepth')
ax3.set_ylabel('Length [m]')
ax3.set_xlabel('Time [days]')
ax3.set_ylim([0,xMax])

cp = ax4.contourf(timeSteps*dt/86400,x,spdLength,
                np.linspace(0,.1,11),
                cmap='cmo.speed',
                extend="both")
cbar = plt.colorbar(cp)
cp2 = ax4.contour(timeSteps*dt/86400,x,spdLength,
                [0.04],
                colors='red',
                linewidths=1,
                linestyles='--',
                alpha=.5)
cbar.set_label('[m/s]')
ax4.grid(alpha=.5)
ax4.set_title('Speed')
ax4.set_ylabel('Length [m]')
ax4.set_xlabel('Time [days]')
ax4.set_ylim([0,xMax])


plt.suptitle(f'Naive depth averaging\nAbove {-zDepth} m')

plt.savefig('figs/momLengthView%s.png' %fileEnding, format='png',dpi=plotDPI)
if(args.silent == 0):
    plt.show()
plt.close()