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
parser.add_argument('-n','--numFrames', nargs=1, default=[250],type=int,
                    help='Max number of samples from time series [defaut = 150]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
args = parser.parse_args()


# Pick cross section to view from file or default
plotDPI = 150
manualMax = None

folders = ['.']
labels = ['glaciome']
# folders = ['hotel_fringe_sgd1300','foxtrot_s_sgd1300']
# labels = ['fringe','none']
# folders = ['hotel_Umin03_sgd1300','foxtrot_s_sgd1300','hotel_Umin05_sgd1300']
# labels = ['0.03','0.04','0.05']
# folders = ['hotel_zhaoMelt_sgd1300','hotel_4x_sgd1300','foxtrot_s_sgd1300']
# labels = ['Zhao','4x','1x']
colors = ['xkcd:blue','xkcd:gray','xkcd:red']
resultFolder = '/results'
fileEnding = ""

# manualMax = 4957200

## Umin melt thresholds 
# thresholds = [0.03,0.04,0.05]
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
    Mr68OverTime = np.zeros(np.shape(timeSteps))
    Mr32OverTime = np.zeros(np.shape(timeSteps))
    for i in range(len(timeSteps)):
        data = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), timeSteps[i])
        fwOverTime[i] = np.nansum(data[0,:,:,:])
        data[2, data[0,:,:,:] == 0] = np.nan #nan all zero melt cells
        MROverTime[i] = np.nanmean(data[2,:,:,:]) #melt ratem m/day
        Mr68OverTime[i] = np.nanpercentile(data[2,:,:,:],68)
        Mr32OverTime[i] = np.nanpercentile(data[2,:,:,:],32)
    ax1.plot(timeSteps*dt/86400,fwOverTime,label=labels[j],color=colors[j])
# Wrap up after plotting everything
ax1.grid(alpha=.5)
# ax1_1=ax1.twinx()
# ax1_1.plot(timeSteps*dt/86400,MROverTime,label=labels[j],color='xkcd:rose',linestyle='-')
# ax1_1.plot(timeSteps*dt/86400,Mr68OverTime,label=labels[j],color='xkcd:light rose',linestyle='-')
# ax1_1.plot(timeSteps*dt/86400,Mr32OverTime,label=labels[j],color='xkcd:light rose',linestyle='-')
ax1_1.set_ylabel('Avg Melt Rate flux [m/d]')
ax1.legend()
ax1.set_title('FW Flux From Bergs Over Time')
ax1.set_ylabel('FW flux [m^3/s]')
ax1.set_xlabel('Time [days]')
# ax1.savefig(str, format='png',dpi=plotDPI)
# plt.show()
# plt.close()

#Open figure before looping
## Plotting water production as function of depth
for i in range(len(folders)):
    folder = folders[i]

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
    
    if(manualMax != None):
        maxStep = manualMax #if manual adjusting needed
    z = np.squeeze(mds.rdmds("%s%s/RC" %(folder,resultFolder))) #midpoints of cells
    dz = np.load("%s/input/dz.npy" %(folder)) # we now save this upon grid generation, much better

    dynName = ['BRGFlx']
    name = ['BRGfwFlx']
    units = ["[m^3/s]"]

    dataStart = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), startStep)
    dataEnd = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), maxStep)
    dataMeltOnly = dataEnd[0,:,:,:] # Only compute spd for melting grids
    data = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'dynDiag'), maxStep)
    spd = (data[2,:,:,:]**2 + data[3,:,:,:]**2 + data[4,:,:,:]**2)**(0.5)
    #above/below threshold
    meltBelowSpd = dataMeltOnly.copy()
    meltBelowSpd[spd > thresholds[i]] = np.nan
    meltAboveSpd = dataMeltOnly.copy()
    meltAboveSpd[spd < thresholds[i]] = np.nan
    # ax2.plot(np.nansum(dataStart[0,:,:,:],axis=(1,2))/dz,z,linestyle='--',color=colors[i],\
    #         label=labels[i] + " 12 hr",alpha=.5)
    ax2.plot(np.nansum(dataEnd[0,:,:,:],axis=(1,2))/dz,z,label=labels[i],color=colors[i])
    # ax2.plot(np.nansum(meltBelowSpd[:,:,:],axis=(1,2))/dz,z,label='spd < $U_{min}$',color=colors[i])
    ax2.plot(np.nansum(meltAboveSpd[:,:,:],axis=(1,2))/dz,z,label='spd > $U_{min}$',color=colors[i],linestyle='--')
# Wrap up after plotting everything

ax2.grid(alpha=.5)
ax2.set_title('Fresh Water Flux Over Depth')
ax2.set_ylim([-450, 0])
ax2.set_ylabel('Depth [m]')
ax2.set_xlabel('Fresh water flux per unit depth [m^2/s]')
ax2.legend()
# plt.savefig(str, format='png',dpi=plotDPI)
# plt.show()
# plt.close()

#Historgram of meltspeeds
ax3_2 = ax3.twinx()
binList = np.arange(0,.20,.20/50)
zDeep = 27 #deepest Drafts of melange

ax3_2.plot([],[],label=labels[0],color=colors[0],linewidth=2)
# ax3_2.plot([],[],label=labels[1],color=colors[1],linewidth=2)

ax3.plot([0,1],[0,0],color='black',linewidth=0.75)
ax3.grid(alpha=.5)
ax3.plot([0.04,0.04],[0,.3],linestyle='-.',color='xkcd:light red',alpha=0.5)

folderData = [None] * (len(folders))
for j in range(len(folders)):
    folder = folders[j]

    x = mds.rdmds("%s%s/XC" %(folder,resultFolder))
    xSlice = np.argmin(np.abs(x[0,:]))

    data = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'dynDiag'), maxStep)
    dataMelt = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'BRGFlx'), maxStep)
    dataMeltOnly = dataMelt[0,:,:,:] # Only compute spd for melting grids
    spd = (data[2,:,:,:]**2 + data[3,:,:,:]**2 + data[4,:,:,:]**2)**(0.5)
    spd[dataMeltOnly == 0] = np.nan
    folderData[j] = spd.flatten()
    localData = spd.flatten()
    wght = np.zeros_like(localData) + 1. / np.count_nonzero(~np.isnan(localData))
    print('length without nan',np.count_nonzero(~np.isnan(localData)))
    print('length with nan', len(localData[localData != np.nan]))
    print(np.sum(wght))
    ax3.hist(localData,bins=binList,label=labels[j],color=colors[j],
        histtype='step',weights=wght,linewidth=2)
    print("50/90th percentile",np.nanmedian(spd),np.nanpercentile(spd,90))
# bigData = [folderData[1],folderData[0],h_data_m,h_data]
bigData = folderData[0]
bplot = ax3_2.boxplot(bigData,vert=False,whis=[10,90],showfliers=False,widths=.75,patch_artist=True)
# bcolors = [colors[1],colors[0],'xkcd:dark grey','xkcd:grey']
bcolors = [colors[0]]
for patch, color in zip(bplot['boxes'], bcolors):
    patch.set_facecolor(color)

ax3.set_ylim([-.06,.230])
ax3.set_title('Histogram of Ocean Speed contributing to melt')
ax3.set_ylabel('Relative abundance []')
ax3.set_xlabel('Speed [m/s]')
ax3.set_yticks([0,.05,.1,.15,.2])
ax3.set_xlim([0,0.20])
ax3_2.set_ylim([.25,22])
ax3_2.set_yticks([])
ax3_2.legend()

#Plot of speed vs depth
# plt.figure()
ax4.plot([0.04,0.04],[10,-450],linestyle='-.',color='xkcd:light red',alpha=0.5)
for j in range(len(folders)):
    folder = folders[j]
    dataMelt = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'BRGFlx'), maxStep)
    dataMeltOnly = dataMelt[0,:,:,:] # Only compute spd for melting grids
    data = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'dynDiag'), maxStep)
    spd = (data[2,:,:,:]**2 + data[3,:,:,:]**2 + data[4,:,:,:]**2)**(0.5)
    spd[dataMeltOnly == 0] = np.nan
    if(j == 0):
        ax4.plot(np.nanmean(spd[:,:,:],axis=(1,2)),z[:],label=labels[j],color=colors[j])
    else:
        ax4.plot(np.nanmean(spd[:,:,:],axis=(1,2)),z[:],label=labels[j],color=colors[j])
str = "figs/spdDepth%s.png" %fileEnding
ax4.grid(alpha=.5)
ax4.set_title('Average Ocean Speed contributing to melt')
ax4.set_ylabel('Depth [m]')
ax4.set_xlabel('Speed [m/s]')
ax4.legend()
ax4.set_ylim([-450, 0])

for label,ax in zip(figLabels,axes.flatten()):
    ax.text(
        ax.get_xlim()[0], ax.get_ylim()[1], label,
        fontsize='x-large', va='bottom',ha='center', fontfamily='sans serif')

plt.savefig('figs/Figure8%s.png' %fileEnding, format='png',dpi=plotDPI)
if(args.silent == 0):
    plt.show()
plt.close()
   