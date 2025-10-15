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
    fileEnding = fileEnding + f'_z{np.abs(zDepth)}'
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
fig, axes = plt.subplots(3, 2, figsize=(12, 8), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]
ax5 = axes[2,0]
ax6 = axes[2,1]
#Total melt over time

for j in range(len(folders)):
    folder = folders[j]

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

    dynName = ['BRGFlx']
    name = ['BRGfwFlx']
    units = ["[m^3/s]"]


    weightsTmp = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), startStep)
    wghts = weightsTmp[0,:,:,:]
    for i in range(len(dz)):
        wghts[i,:,:] = dz[i]

    timeSteps = np.arange(startStep, maxStep + 1, sizeStep)
    fwOverTime = np.zeros(np.shape(timeSteps))
    MROverTime = np.zeros(np.shape(timeSteps))
    Mr25Time = np.zeros(np.shape(timeSteps))
    Mr32Time = np.zeros(np.shape(timeSteps))
    Mr50Time = np.zeros(np.shape(timeSteps))
    Mr68Time = np.zeros(np.shape(timeSteps))
    Mr75Time = np.zeros(np.shape(timeSteps))
    Mr95Time = np.zeros(np.shape(timeSteps))
    TavgTime = np.zeros(np.shape(timeSteps))
    T25Time = np.zeros(np.shape(timeSteps))
    T32Time = np.zeros(np.shape(timeSteps))
    T50Time = np.zeros(np.shape(timeSteps))
    T68Time = np.zeros(np.shape(timeSteps))
    T75Time = np.zeros(np.shape(timeSteps))
    T95Time = np.zeros(np.shape(timeSteps))
    UavgTime = np.zeros(np.shape(timeSteps))
    U25Time = np.zeros(np.shape(timeSteps))
    U32Time = np.zeros(np.shape(timeSteps))
    U50Time = np.zeros(np.shape(timeSteps))
    U68Time = np.zeros(np.shape(timeSteps))
    U75Time = np.zeros(np.shape(timeSteps))
    U95Time = np.zeros(np.shape(timeSteps))
    SavgTime = np.zeros(np.shape(timeSteps))
    S25Time = np.zeros(np.shape(timeSteps))
    S32Time = np.zeros(np.shape(timeSteps))
    S50Time = np.zeros(np.shape(timeSteps))
    S68Time = np.zeros(np.shape(timeSteps))
    S75Time = np.zeros(np.shape(timeSteps))
    S95Time = np.zeros(np.shape(timeSteps))
    WavgTime = np.zeros(np.shape(timeSteps))
    W25Time = np.zeros(np.shape(timeSteps))
    W32Time = np.zeros(np.shape(timeSteps))
    W50Time = np.zeros(np.shape(timeSteps))
    W68Time = np.zeros(np.shape(timeSteps))
    W75Time = np.zeros(np.shape(timeSteps))
    W95Time = np.zeros(np.shape(timeSteps))
    # Depth view
    fwDepth = np.zeros([len(z),len(timeSteps)])
    bDepth = np.zeros([len(z),len(timeSteps)])
    tDepth = np.zeros([len(z),len(timeSteps)])
    sDepth = np.zeros([len(z),len(timeSteps)])
    uDepth = np.zeros([len(z),len(timeSteps)])
    wDepth = np.zeros([len(z),len(timeSteps)])
    pDepth = np.zeros([len(z),len(timeSteps)])
    hfDepth = np.zeros([len(z),len(timeSteps)])
    # Length view
    fwLength = np.zeros([len(x),len(timeSteps)])
    bLength = np.zeros([len(x),len(timeSteps)])
    tLength = np.zeros([len(x),len(timeSteps)])
    sLength = np.zeros([len(x),len(timeSteps)])
    uLength = np.zeros([len(x),len(timeSteps)])
    wLength = np.zeros([len(x),len(timeSteps)])
    pLength = np.zeros([len(x),len(timeSteps)])
    hfLength = np.zeros([len(x),len(timeSteps)])
    for i in range(len(timeSteps)):
        data = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), timeSteps[i])
        dataOcean = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'dynDiag'), timeSteps[i])
        dataPress = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'presDiag'), timeSteps[i])
        fwOverTime[i] = np.nansum(data[0,:,:,:])
        melangeMask = data[0,:,:,:].copy()
        melangeMask[melangeMask != 0] = 1
        melangeMask[data[0,:,:,:] < -1] = 1 #attemp to recapture freezing limit ice
        melangeMask[:,np.nansum(data[0,:,:,:],axis=0) == 0] = 0 #attemp to recapture freezing limit ice
        melangeMask[:,:,1] = 0 #exclude glacier face and plume
        data[:, melangeMask == 0] = np.nan #nan all zero melt cells
        dataOcean[:, melangeMask == 0] = np.nan #nan all zero melt cells# spd = ((dataOcean[2,:,:,:]**2 + dataOcean[3,:,:,:]**2 + dataOcean[4,:,:,:]**2)**(.5))
        meanP = np.nanmean(dataPress[0,:,int(ny/2),2:-2],axis=(1)) #find anomoly relative to centerline average
        # print(meanP)
        dataPress[0,:,:,:] = dataPress[0,:,:,:] - meanP[:,None,None] #get anomoly before NANing
        dataPress[:, melangeMask == 0] = np.nan
        spd = dataOcean[2,:,:,:]
        if(sgd):
            sgdFactor = f(i*dt)/1000
        else:
            sgdFactor = 1
        #depth values
        fwDepth[:,i] = np.nansum(data[0,:,:,xSlice:],axis=(1,2))/dz/sgdFactor
        bDepth[:,i] = np.nanmean(data[2,:,:,xSlice:],axis=(1,2))/sgdFactor
        tDepth[:,i] = np.nanmean(dataOcean[0,:,:,xSlice:],axis=(1,2))
        sDepth[:,i] = np.nanmean(dataOcean[1,:,:,xSlice:],axis=(1,2))
        uDepth[:,i] = np.nanmean(spd[:,:,xSlice:],axis=(1,2))/sgdFactor
        wDepth[:,i] = np.nanmean(dataOcean[3,:,:,xSlice:],axis=(1,2))/sgdFactor
        pDepth[:,i] = np.nanmean(dataPress[0,:,:,xSlice:],axis=(1,2))/sgdFactor
        hfDepth[:,i] = np.nanmean(data[5,:,:,xSlice:],axis=(1,2))
        #length values
        
        fwLength[:,i] = np.nansum(data[0,:zSlice,:,:],axis=(0,1))/(x[1]-x[0])/sgdFactor
        bLength[:,i] = np.nanmean(data[2,:zSlice,:,:],axis=(0,1))/sgdFactor
        tLength[:,i] = np.nanmean(dataOcean[0,:zSlice,:,:],axis=(0,1))
        sLength[:,i] = np.nanmean(dataOcean[1,:zSlice,:,:],axis=(0,1))
        uLength[:,i] = np.nanmean(spd[:zSlice,:,:],axis=(0,1))/sgdFactor
        wLength[:,i] = np.nanmean(dataOcean[3,:zSlice,:,:],axis=(0,1))/sgdFactor
        pLength[:,i] = np.nanmean(dataPress[0,:zSlice,:,:],axis=(0,1))/sgdFactor
        hfLength[:,i] = np.nanmean(data[5,:zSlice,:,:],axis=(0,1))
        #percentile values
        MROverTime[i] = np.nanmean(data[2,:,:,:]) #melt ratem m/day
        Mr25Time[i] = np.nanpercentile(data[2,:,:,:],25,weights=wghts,method='inverted_cdf')
        Mr32Time[i] = np.nanpercentile(data[2,:,:,:],32,weights=wghts,method='inverted_cdf')
        Mr50Time[i] = np.nanpercentile(data[2,:,:,:],50,weights=wghts,method='inverted_cdf')
        Mr68Time[i] = np.nanpercentile(data[2,:,:,:],68,weights=wghts,method='inverted_cdf')
        Mr75Time[i] = np.nanpercentile(data[2,:,:,:],75,weights=wghts,method='inverted_cdf')
        Mr95Time[i] = np.nanpercentile(data[2,:,:,:],95,weights=wghts,method='inverted_cdf')
        TavgTime[i] = np.nanmean(dataOcean[0,:,:,:]) #temp 
        T25Time[i] = np.nanpercentile(dataOcean[0,:,:,:],25,weights=wghts,method='inverted_cdf')
        T32Time[i] = np.nanpercentile(dataOcean[0,:,:,:],32,weights=wghts,method='inverted_cdf')
        T50Time[i] = np.nanpercentile(dataOcean[0,:,:,:],50,weights=wghts,method='inverted_cdf')
        T68Time[i] = np.nanpercentile(dataOcean[0,:,:,:],68,weights=wghts,method='inverted_cdf')
        T75Time[i] = np.nanpercentile(dataOcean[0,:,:,:],75,weights=wghts,method='inverted_cdf')
        T95Time[i] = np.nanpercentile(dataOcean[0,:,:,:],95,weights=wghts,method='inverted_cdf')
        SavgTime[i] = np.nanmean(dataOcean[1,:,:,:]) #salt 
        S25Time[i] = np.nanpercentile(dataOcean[1,:,:,:],25,weights=wghts,method='inverted_cdf')
        S32Time[i] = np.nanpercentile(dataOcean[1,:,:,:],32,weights=wghts,method='inverted_cdf')
        S50Time[i] = np.nanpercentile(dataOcean[1,:,:,:],50,weights=wghts,method='inverted_cdf')
        S68Time[i] = np.nanpercentile(dataOcean[1,:,:,:],68,weights=wghts,method='inverted_cdf')
        S75Time[i] = np.nanpercentile(dataOcean[1,:,:,:],75,weights=wghts,method='inverted_cdf')
        S95Time[i] = np.nanpercentile(dataOcean[1,:,:,:],95,weights=wghts,method='inverted_cdf')
        WavgTime[i] = np.nanmean(dataOcean[3,:,:,:]) #vert speed 
        W25Time[i] = np.nanpercentile(dataOcean[3,:,:,:],25,weights=wghts,method='inverted_cdf')
        W32Time[i] = np.nanpercentile(dataOcean[3,:,:,:],32,weights=wghts,method='inverted_cdf')
        W50Time[i] = np.nanpercentile(dataOcean[3,:,:,:],50,weights=wghts,method='inverted_cdf')
        W68Time[i] = np.nanpercentile(dataOcean[3,:,:,:],68,weights=wghts,method='inverted_cdf')
        W75Time[i] = np.nanpercentile(dataOcean[3,:,:,:],75,weights=wghts,method='inverted_cdf')
        W95Time[i] = np.nanpercentile(dataOcean[3,:,:,:],95,weights=wghts,method='inverted_cdf')
        UavgTime[i] = np.nanmean(spd) #speed 
        U25Time[i] = np.nanpercentile(spd,25,weights=wghts,method='inverted_cdf')
        U32Time[i] = np.nanpercentile(spd,32,weights=wghts,method='inverted_cdf')
        U50Time[i] = np.nanpercentile(spd,50,weights=wghts,method='inverted_cdf')
        U68Time[i] = np.nanpercentile(spd,68,weights=wghts,method='inverted_cdf')
        U75Time[i] = np.nanpercentile(spd,75,weights=wghts,method='inverted_cdf')
        U95Time[i] = np.nanpercentile(spd,95,weights=wghts,method='inverted_cdf')
    ax1.plot(timeSteps*dt/86400,fwOverTime,label=labels[j],color=colors[j])
    # Wrap up after plotting everything
    ax1.grid(alpha=.5)
    ax1.legend()
    ax1.set_title('FW Flux From Bergs Over Time')
    ax1.set_ylabel('FW flux [m^3/s]')
    ax1.set_xlabel('Time [days]')
    # ax1.savefig(str, format='png',dpi=plotDPI)
    # plt.show()
    # plt.close()


    ax2.plot(timeSteps*dt/86400,Mr95Time,color='xkcd:light purple',linestyle='--',label='95%')
    ax2.plot(timeSteps*dt/86400,Mr25Time,color='xkcd:light purple',linestyle='-',label='$\\pm$25%')
    ax2.plot(timeSteps*dt/86400,Mr32Time,color='xkcd:purple',linestyle='--',label='$\\pm$18%')
    ax2.plot(timeSteps*dt/86400,Mr50Time,color='xkcd:dark purple',linestyle='-',label='Median')
    ax2.plot(timeSteps*dt/86400,MROverTime,color='xkcd:gray',linestyle='-',label='Average')
    ax2.plot(timeSteps*dt/86400,Mr68Time,color='xkcd:purple',linestyle='--')
    ax2.plot(timeSteps*dt/86400,Mr75Time,color='xkcd:light purple',linestyle='-')

    ax2.grid(alpha=.5)
    ax2.legend()
    ax2.set_title('Melt Rates Over Time')
    ax2.set_ylabel('Melt Rate flux [m/d]')
    ax2.set_xlabel('Time [days]')

    ax3.plot(timeSteps*dt/86400,T95Time,color='xkcd:light orange',linestyle='--',label='95%')
    ax3.plot(timeSteps*dt/86400,T25Time,color='xkcd:light orange',linestyle='-',label='$\\pm$25%')
    ax3.plot(timeSteps*dt/86400,T32Time,color='xkcd:orange',linestyle='--',label='$\\pm$18%')
    ax3.plot(timeSteps*dt/86400,T50Time,color='xkcd:dark orange',linestyle='-',label='Median')
    ax3.plot(timeSteps*dt/86400,TavgTime,color='xkcd:gray',linestyle='-',label='Average')
    ax3.plot(timeSteps*dt/86400,T68Time,color='xkcd:orange',linestyle='--')
    ax3.plot(timeSteps*dt/86400,T75Time,color='xkcd:light orange',linestyle='-')

    ax3.grid(alpha=.5)
    ax3.legend()
    ax3.set_title('Temperature Over Time')
    ax3.set_ylabel('Temp [C]')
    ax3.set_xlabel('Time [days]')

    ax4.plot(timeSteps*dt/86400,U95Time,color='xkcd:light blue',linestyle='--',label='95%')
    ax4.plot(timeSteps*dt/86400,U25Time,color='xkcd:light blue',linestyle='-',label='$\\pm$25%')
    ax4.plot(timeSteps*dt/86400,U32Time,color='xkcd:blue',linestyle='--',label='$\\pm$18%')
    ax4.plot(timeSteps*dt/86400,U50Time,color='xkcd:dark blue',linestyle='-',label='Median')
    ax4.plot(timeSteps*dt/86400,UavgTime,color='xkcd:gray',linestyle='-',label='Average')
    ax4.plot(timeSteps*dt/86400,U68Time,color='xkcd:blue',linestyle='--')
    ax4.plot(timeSteps*dt/86400,U75Time,color='xkcd:light blue',linestyle='-')

    ax4.grid(alpha=.5)
    ax4.legend()
    ax4.set_title('Ocean Speed Over Time')
    ax4.set_ylabel('Speed [m/s]')
    ax4.set_xlabel('Time [days]')

    ax5.plot(timeSteps*dt/86400,S95Time,color='xkcd:light green',linestyle='--',label='95%')
    ax5.plot(timeSteps*dt/86400,S25Time,color='xkcd:light green',linestyle='-',label='$\\pm$25%')
    ax5.plot(timeSteps*dt/86400,S32Time,color='xkcd:green',linestyle='--',label='$\\pm$18%')
    ax5.plot(timeSteps*dt/86400,S50Time,color='xkcd:dark green',linestyle='-',label='Median')
    ax5.plot(timeSteps*dt/86400,SavgTime,color='xkcd:gray',linestyle='-',label='Average')
    ax5.plot(timeSteps*dt/86400,S68Time,color='xkcd:green',linestyle='--')
    ax5.plot(timeSteps*dt/86400,S75Time,color='xkcd:light green',linestyle='-')

    ax5.grid(alpha=.5)
    ax5.legend()
    ax5.set_title('Ocean Salinity Over Time')
    ax5.set_ylabel('Salt [PSU]')
    ax5.set_xlabel('Time [days]')

    ax6.plot(timeSteps*dt/86400,W95Time,color='xkcd:light red',linestyle='--',label='95%')
    ax6.plot(timeSteps*dt/86400,W25Time,color='xkcd:light red',linestyle='-',label='$\\pm$25%')
    ax6.plot(timeSteps*dt/86400,W32Time,color='xkcd:red',linestyle='--',label='$\\pm$18%')
    ax6.plot(timeSteps*dt/86400,W50Time,color='xkcd:dark red',linestyle='-',label='Median')
    ax6.plot(timeSteps*dt/86400,WavgTime,color='xkcd:gray',linestyle='-',label='Average')
    ax6.plot(timeSteps*dt/86400,W68Time,color='xkcd:red',linestyle='--')
    ax6.plot(timeSteps*dt/86400,W75Time,color='xkcd:light red',linestyle='-')

    ax6.grid(alpha=.5)
    ax6.legend()
    ax6.set_title('Ocean W Over Time')
    ax6.set_ylabel('W [m/s]')
    ax6.set_xlabel('Time [days]')

plt.savefig('figs/meltDetails%s.png' %fileEnding, format='png',dpi=plotDPI)
if(args.silent == 0):
    plt.show()
plt.close()

## Next figure
fig, axes = plt.subplots(4, 2, figsize=(12, 8), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]
ax5 = axes[2,0]
ax6 = axes[2,1]
ax7 = axes[3,0]
ax8 = axes[3,1]

iceRange = np.linspace(0,1,11)
pRange = np.linspace(-.4,.4,31)

#minor cleaning, nansum returns 0 if all nan, nanmean returns nan if all nan. We prefer the 2nd behavior
fwDepth[np.isnan(tDepth)] = np.nan

zMin=np.min(z[~np.isnan(np.nanmean(fwDepth,axis=1))])
cp = ax1.contourf(timeSteps*dt/86400,z,fwDepth,
                cmap='cmo.dense')
cbar = plt.colorbar(cp)
cbar.set_label('FW flux/meter [m^2/s]')
ax1.grid(alpha=.5)
ax1.set_title('FW Flux (un-normalized SUM)')
ax1.set_ylabel('Depth [m]')
ax1.set_xlabel('Time [days]')
ax1.set_ylim([zMin,0])

cp = ax2.contourf(timeSteps*dt/86400,z,bDepth,
                np.linspace(0,2,31),
                cmap=meltCmap,
                extend='both')
cbar = plt.colorbar(cp)
cbar.set_label('Averge Meltrate [m/day]')
ax2.grid(alpha=.5)
ax2.set_title('Meltrate')
ax2.set_ylabel('Depth [m]')
ax2.set_xlabel('Time [days]')
ax2.set_ylim([zMin,0])

cp = ax3.contourf(timeSteps*dt/86400,z,tDepth,
                tempRange,
                cmap=tempCmap)
cbar = plt.colorbar(cp)
cbar.set_label('Temperature [C]')
ax3.grid(alpha=.5)
ax3.set_title('Temperature')
ax3.set_ylabel('Depth [m]')
ax3.set_xlabel('Time [days]')
ax3.set_ylim([zMin,0])

cp = ax4.contourf(timeSteps*dt/86400,z,uDepth,
                uRange,
                cmap=uCmap)
cbar = plt.colorbar(cp)
cbar.set_label('Velocity [m/s]')
ax4.grid(alpha=.5)
ax4.set_title('U Vel Down Fjord')
ax4.set_ylabel('Depth [m]')
ax4.set_xlabel('Time [days]')
ax4.set_ylim([zMin,0])

cp = ax5.contourf(timeSteps*dt/86400,z,sDepth,
                saltRange,
                cmap=saltCmap)
cbar = plt.colorbar(cp)
cbar.set_label('Salinity [PSU]')
ax5.grid(alpha=.5)
ax5.set_title('Salinity')
ax5.set_ylabel('Depth [m]')
ax5.set_xlabel('Time [days]')
ax5.set_ylim([zMin,0])

cp = ax6.contourf(timeSteps*dt/86400,z,wDepth,
                wRange,
                cmap=wCmap,
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('Upwelling [m/s]')
ax6.grid(alpha=.5)
ax6.set_title('Upwelling')
ax6.set_ylabel('Depth [m]')
ax6.set_xlabel('Time [days]')
ax6.set_ylim([zMin,0])

cp = ax7.contourf(timeSteps*dt/86400,z,pDepth,
                pRange,
                cmap="PuOr_r",
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('Pressure [m^2/s^2]')
ax7.grid(alpha=.5)
ax7.set_title('HYD Pressure Anom')
ax7.set_ylabel('Depth [m]')
ax7.set_xlabel('Time [days]')
ax7.set_ylim([zMin,0])

cp = ax8.contourf(timeSteps*dt/86400,z,1-hfDepth,
                iceRange,
                cmap="Grays",
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('Ice Vol Fraction [ ]')
ax8.grid(alpha=.5)
ax8.set_title('Ice Fill')
ax8.set_ylabel('Depth [m]')
ax8.set_xlabel('Time [days]')
ax8.set_ylim([zMin,0])

plt.suptitle(f'View Beyond x = {xReach} m')
plt.savefig('figs/meltDepthView%s.png' %fileEnding, format='png',dpi=plotDPI)
if(args.silent == 0):
    plt.show()
plt.close()

## Next figure
fig, axes = plt.subplots(4, 2, figsize=(12, 8), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]
ax5 = axes[2,0]
ax6 = axes[2,1]
ax7 = axes[3,0]
ax8 = axes[3,1]

#minor cleaning, nansum returns 0 if all nan, nanmean returns nan if all nan. We prefer the 2nd behavior
fwLength[np.isnan(tLength)] = np.nan

xMax=np.max(x[~np.isnan(np.nanmean(fwLength,axis=1))])
cp = ax1.contourf(timeSteps*dt/86400,x,fwLength,
                cmap='cmo.dense')
cbar = plt.colorbar(cp)
cbar.set_label('FW flux/meter [m^2/s]')
ax1.grid(alpha=.5)
ax1.set_title('FW Flux (un-normalized SUM)')
ax1.set_ylabel('Length [m]')
ax1.set_xlabel('Time [days]')
ax1.set_ylim([0,xMax])

cp = ax2.contourf(timeSteps*dt/86400,x,bLength,
                np.linspace(0,2,31),
                cmap=meltCmap,
                extend='both')
cbar = plt.colorbar(cp)
cbar.set_label('Averge Meltrate [m/day]')
ax2.grid(alpha=.5)
ax2.set_title('Meltrate')
ax2.set_ylabel('Length [m]')
ax2.set_xlabel('Time [days]')
ax2.set_ylim([0,xMax])

cp = ax3.contourf(timeSteps*dt/86400,x,tLength,
                tempRange,
                cmap=tempCmap)
cbar = plt.colorbar(cp)
cbar.set_label('Temperature [C]')
ax3.grid(alpha=.5)
ax3.set_title('Temperature')
ax3.set_ylabel('Length [m]')
ax3.set_xlabel('Time [days]')
ax3.set_ylim([0,xMax])

cp = ax4.contourf(timeSteps*dt/86400,x,uLength,
                uRange,
                cmap=uCmap)
cbar = plt.colorbar(cp)
cbar.set_label('Velocity [m/s]')
ax4.grid(alpha=.5)
ax4.set_title('U Vel Down Fjord')
ax4.set_ylabel('Length [m]')
ax4.set_xlabel('Time [days]')
ax4.set_ylim([0,xMax])

cp = ax5.contourf(timeSteps*dt/86400,x,sLength,
                saltRange,
                cmap=saltCmap)
cbar = plt.colorbar(cp)
cbar.set_label('Salinity [PSU]')
ax5.grid(alpha=.5)
ax5.set_title('Salinity')
ax5.set_ylabel('Length [m]')
ax5.set_xlabel('Time [days]')
ax5.set_ylim([0,xMax])

cp = ax6.contourf(timeSteps*dt/86400,x,wLength,
                wRange,
                cmap=wCmap,
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('Upwelling [m/s]')
ax6.grid(alpha=.5)
ax6.set_title('Upwelling')
ax6.set_ylabel('Length [m]')
ax6.set_xlabel('Time [days]')
ax6.set_ylim([0,xMax])

cp = ax7.contourf(timeSteps*dt/86400,x,pLength,
                pRange,
                cmap="PuOr_r",
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('Pressure [m^2/s^2]')
ax7.grid(alpha=.5)
ax7.set_title('HYD Pressure Anom')
ax7.set_ylabel('Length [m]')
ax7.set_xlabel('Time [days]')
ax7.set_ylim([0,xMax])

cp = ax8.contourf(timeSteps*dt/86400,x,1-hfLength,
                iceRange,
                cmap="Grays",
                extend="both")
cbar = plt.colorbar(cp)
cbar.set_label('Ice Fraction [ ]')
ax8.grid(alpha=.5)
ax8.set_title('Ice Frac')
ax8.set_ylabel('Length [m]')
ax8.set_xlabel('Time [days]')
ax8.set_ylim([0,xMax])

plt.suptitle(f'Naive depth averaging\nAbove {-zDepth} m')

plt.savefig('figs/meltLengthView%s.png' %fileEnding, format='png',dpi=plotDPI)
if(args.silent == 0):
    plt.show()
plt.close()