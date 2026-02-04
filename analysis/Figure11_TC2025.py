from MITgcmutils import mds
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import xarray as xr
import fileinput
from matplotlib.colors import LinearSegmentedColormap

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 20000
zDepth = -10
plotDPI = 300
cleanPNGs = True

def red_yellow_white_cyan_blue():
    N = 256
    cols = ['#000055', '#0000f5', '#008cff',
            # '#7affff', '#ffffff', '#ffff83',
            '#9ee6e6', '#ececec', '#eaea99',
            '#ff9d0c', '#ff0500', '#5f0000']
    return LinearSegmentedColormap.from_list('custom', cols, N)


folders = ['FJ_NoIcebergs','FJ_NoBlockDrag','FJ']
titles = ['FJ_NoIcebergs', 'FJ_NoBlockDrag', 'FJ']
colors = ['xkcd:pumpkin','xkcd:green','xkcd:blue']
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


zSlice = 15
stepsToAverage = 10

# fig = plt.figure(figsize=(12, 6))
# ax1 = plt.subplot(2,3,1)
# ax2 = plt.subplot(2,3,2)
# ax3 = plt.subplot(2,3,3)
# ax4 = plt.subplot(2,3,4)
# ax5 = plt.subplot(2,3,5)
# ax6 = plt.subplot(2,3,6)
# axes = [ax1, ax2, ax3, ax4, ax5, ax6]
labels = ["(a)","(b)","(c)","(d)","(e)","(f)"]
# labels = ["E","F","G"]
fig, axes = plt.subplots(3, 2, figsize=(12, 4), layout="constrained",width_ratios=[10, 1])
ax1 = axes[0,0]
ax2 = axes[1,0]
ax3 = axes[2,0]
ax4 = axes[0,1]
ax5 = axes[1,1]
ax6 = axes[2,1]
# fig, axes1 = plt.subplots(3, 1, figsize=(2.5, 5.25), layout="constrained")
# axes=np.array([[None, axes1[0]],[None, axes1[1]],[None,axes1[2]]])
# print(axes)

lvl = np.linspace(-.10, .10, 127)
# plt.subplot(131)

#Find time steps to take        
maxStep = 0
sizeStep = 1e10
startStep = 1e10
frameLimit = 100
for i in range(len(folders)):
    folder = folders[i]
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
            if((maxStep-startStep)/sizeStep > frameLimit):   #if more than # frames, downscale to be less than #
                dwnScale = int(np.ceil(((maxStep-startStep)/sizeStep)/frameLimit))
                print('Reducing time resolution by', dwnScale)
                sizeStep = sizeStep * dwnScale

    print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

    dt = 0.0   
    for line in fileinput.input('%s/input/data' %folders[0]):
            if "deltaT=" in line:
                dt = float(line[8:-2])
    print('dt is loaded as', dt)

    # Plume
    y = np.squeeze(mds.rdmds("%s%s/YC"%(folder,resultFolder))[:,0]-3400)
    x = np.squeeze(mds.rdmds("%s%s/XC"%(folder,resultFolder))[0,:])
    z = np.squeeze(mds.rdmds("%s%s/RC"%(folder,resultFolder)))

    xSlice = np.argmin(np.abs(x - xCrossSection))
    print('xslice at location %.3f' %x[xSlice])
    print('zslice at depth %.3f' %z[zSlice])


    tmpSize = np.array(np.shape(mds.rdmds("%s%s/dynDiag"%(folder,resultFolder), maxStep)))
    tmpSize[0] = stepsToAverage
    dataArray = np.zeros(tmpSize)
    timeSteps = np.arange(maxStep-(stepsToAverage-1)*sizeStep, maxStep + 1, sizeStep)
    print("averaging over final",(timeSteps[-1]-timeSteps[0])*dt/86400,"days")
    for j in range(len(timeSteps)):
        dataArray[j,:,:,:] = mds.rdmds("%s%s/dynDiag"%(folder,resultFolder), timeSteps[j])[2,:,:,:]
    data = np.nanmean(dataArray,axis=0)
    cp = axes[i,0].pcolormesh(
                x[:]/1e3,
                y[1:-1],
                data[zSlice,1:-1,:],
                vmin=np.min(lvl),
                vmax=np.max(lvl),
                #extend="both",
                cmap=red_yellow_white_cyan_blue(),
                )
    mask = np.zeros_like(data[0,1:-1,:])
    mask[:,xSlice:(2*xSlice)] = 1
    mask[mask==0] = np.nan
    cmask = axes[i,0].pcolormesh(
                x[:]/1e3,
                y[1:-1],
                mask,
                vmin=0,
                vmax=1,
                #extend="both",
                cmap='Grays',
                alpha=0.15
                )
    czero = axes[i,0].contour(
                x[:]/1e3,
                y[1:-1],
                data[zSlice,1:-1,:],
                [0],
                colors='gray',
                linewidths=0.5,
                alpha=0.5)
    axes[i,0].clabel(czero, inline=3, fontsize=8)
    axes[i,0].set_title(titles[i])
    axes[i,0].set_aspect(1/1000)
    if(i == 1):
        axes[i,0].set_ylabel('Across fjord [m]')

    if(i==2):
        axes[i,0].set_xlabel('Along fjord [km]')
        cbar = plt.colorbar(cp,orientation="horizontal",aspect=40)
        cbar.set_label('Along fjord velocity $u$ [m/s]')
    else:
        axes[i,0].set_xticks([])

    cp = axes[i,1].pcolormesh(
            y[1:-1],
            z[:],
            np.nanmean(data[:,1:-1,xSlice:(2*xSlice)],axis=2),
            vmin=np.min(lvl),
            vmax=np.max(lvl),
            #extend="both",
            cmap=red_yellow_white_cyan_blue(),
            )
    axes[i,1].plot([y[1]-200,y[-1]-200],[z[zSlice],z[zSlice]],color='gray',linestyle='--',alpha=0.5)
    czero = axes[i,1].contour(
                y[1:-1],
                z[:],
                np.nanmean(data[:,1:-1,xSlice:(2*xSlice)],axis=2),
                [0],
                colors='gray',
                linewidths=0.5,
                alpha=0.5)
    axes[i,1].clabel(czero, inline=3, fontsize=8)
    axes[i,1].yaxis.tick_right()
    # axes[i,1].set_title(titles[i])
    if(i == 1):
        axes[i,1].set_ylabel('Depth [m]')
    axes[i,1].yaxis.set_label_position("right")
    if(i==2):
        axes[i,1].set_xlabel('Across fjord [m]')
        # cbar = plt.colorbar(cp,orientation="horizontal",aspect=10)
        # cbar.set_label('Along fjord velocity $u$ [m/s]')
    else:
        axes[i,1].set_xticks([])
    if(i==0):
        # axes[i,1].set_title('Cross Section')
        pass
# plt.tight_layout()
for label,ax in zip(labels,axes.reshape(6,1)):
    ax=ax[0]
    ax.text(
        ax.get_xlim()[0], ax.get_ylim()[1], label,
        fontsize='x-large', va='bottom',ha='center', fontfamily='sans serif')

plt.savefig('figs/Figure11%s.png' %fileEnding, format='png',dpi=plotDPI)
plt.show()
plt.close()

