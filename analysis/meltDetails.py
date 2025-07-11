from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import xarray as xr
import argparse

# Take input options. Some defaults are set here, so be aware
parser = argparse.ArgumentParser(description='Plot options for fresh water flux over time')
parser.add_argument('-s','--silent', action='count', default=0,
                    help='Option to silence showing of plots')
parser.add_argument('-n','--numFrames', nargs=1, default=[150],type=int,
                    help='Max number of samples from time series [defaut = 150]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
args = parser.parse_args()


# Pick cross section to view from file or default
plotDPI = 150
manualMax = None

folders = ['.']
labels = ['glaciome']
colors = ['xkcd:blue','xkcd:gray','xkcd:red']
resultFolder = '/results'
fileEnding = ""

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


figLabels = ["(a)","(b)","(c)","(d)","(d)","(f)"]
fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]

#Total melt over time

for j in range(len(folders)):
    folder = folders[j]

    dt = 0.0   
    for line in fileinput.input('%s/input/data' %folder):
            if "deltaT=" in line:
                dt = float(line[8:-2])
    print('dt is loaded as', dt)

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
    for i in range(len(timeSteps)):
        data = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), timeSteps[i])
        dataOcean = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'dynDiag'), timeSteps[i])
        fwOverTime[i] = np.nansum(data[0,:,:,:])
        data[2, data[0,:,:,:] == 0] = np.nan #nan all zero melt cells
        dataOcean[:, data[0,:,:,:] == 0] = np.nan #nan all zero melt cells
        spd = ((dataOcean[2,:,:,:]**2 + dataOcean[3,:,:,:]**2 + dataOcean[4,:,:,:]**2)**(.5))
        MROverTime[i] = np.nanmean(data[2,:,:,:]) #melt ratem m/day
        Mr25Time[i] = np.nanpercentile(data[2,:,:,:],25)
        Mr32Time[i] = np.nanpercentile(data[2,:,:,:],32)
        Mr50Time[i] = np.nanpercentile(data[2,:,:,:],50)
        Mr68Time[i] = np.nanpercentile(data[2,:,:,:],68)
        Mr75Time[i] = np.nanpercentile(data[2,:,:,:],75)
        Mr95Time[i] = np.nanpercentile(data[2,:,:,:],95)
        TavgTime[i] = np.nanmean(dataOcean[0,:,:,:]) #temp 
        T25Time[i] = np.nanpercentile(dataOcean[0,:,:,:],25)
        T32Time[i] = np.nanpercentile(dataOcean[0,:,:,:],32)
        T50Time[i] = np.nanpercentile(dataOcean[0,:,:,:],50)
        T68Time[i] = np.nanpercentile(dataOcean[0,:,:,:],68)
        T75Time[i] = np.nanpercentile(dataOcean[0,:,:,:],75)
        T95Time[i] = np.nanpercentile(dataOcean[0,:,:,:],95)
        UavgTime[i] = np.nanmean(spd) #speed 
        U25Time[i] = np.nanpercentile(spd,25)
        U32Time[i] = np.nanpercentile(spd,32)
        U50Time[i] = np.nanpercentile(spd,50)
        U68Time[i] = np.nanpercentile(spd,68)
        U75Time[i] = np.nanpercentile(spd,75)
        U95Time[i] = np.nanpercentile(spd,95)
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


    ax2.plot(timeSteps*dt/86400,Mr95Time,color='xkcd:light green',linestyle='--',label='95%')
    ax2.plot(timeSteps*dt/86400,Mr25Time,color='xkcd:light green',linestyle='-',label='$\\pm$25%')
    ax2.plot(timeSteps*dt/86400,Mr32Time,color='xkcd:green',linestyle='--',label='$\\pm$18%')
    ax2.plot(timeSteps*dt/86400,Mr50Time,color='xkcd:dark green',linestyle='-',label='Median')
    ax2.plot(timeSteps*dt/86400,MROverTime,color='xkcd:gray',linestyle='-',label='Average')
    ax2.plot(timeSteps*dt/86400,Mr68Time,color='xkcd:green',linestyle='--')
    ax2.plot(timeSteps*dt/86400,Mr75Time,color='xkcd:light green',linestyle='-')

    ax2.grid(alpha=.5)
    ax2.legend()
    ax2.set_title('Melt Rates Over Time')
    ax2.set_ylabel('Melt Rate flux [m/d]')
    ax2.set_xlabel('Time [days]')

    ax3.plot(timeSteps*dt/86400,T95Time,color='xkcd:light red',linestyle='--',label='95%')
    ax3.plot(timeSteps*dt/86400,T25Time,color='xkcd:light red',linestyle='-',label='$\\pm$25%')
    ax3.plot(timeSteps*dt/86400,T32Time,color='xkcd:red',linestyle='--',label='$\\pm$18%')
    ax3.plot(timeSteps*dt/86400,T50Time,color='xkcd:dark red',linestyle='-',label='Median')
    ax3.plot(timeSteps*dt/86400,TavgTime,color='xkcd:gray',linestyle='-',label='Average')
    ax3.plot(timeSteps*dt/86400,T68Time,color='xkcd:red',linestyle='--')
    ax3.plot(timeSteps*dt/86400,T75Time,color='xkcd:light red',linestyle='-')

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

plt.savefig('figs/meltDetails%s.png' %fileEnding, format='png',dpi=plotDPI)
if(args.silent == 0):
    plt.show()
plt.close()
   