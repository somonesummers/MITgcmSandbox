from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import argparse

# Take input options. Some defaults are set here, so be aware
parser = argparse.ArgumentParser(description='Plot options for fjord conditions')
parser.add_argument('-f','--files', nargs='*', default=['.'],
                    help='folders to plot [default is .]')
parser.add_argument('-l','--labels', nargs='*', default=None,
                    help='file string label plots [default = None]')
parser.add_argument('-s','--silent', action='count', default=0,
                    help='Option to silence showing of plots')
parser.add_argument('-n','--numFrames', nargs=1, default=[250],type=int,
                    help='Max number of samples from time series [defaut = 250]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = [500, 860],
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('-dpi','--dpi', nargs='?', default=400,type=int,
                    help='dpi to print [default = 400]')
parser.add_argument('-u','--uRange', nargs='?', default=None,type=float,
                    help='Max U velocity, default to uRange from helperfile')
parser.add_argument('-p','--profile', action='count', default=0,
                    help='option to plot with observed profiles')
args = parser.parse_args()

plt.rcParams.update({'font.size': 14})

# Should start in march, end 1 year later

t_range = [-2.2, 4.2]
if(args.uRange == None):
    u_range = [-0.1, 0.1]
else:
    u_range = [-args.uRange, args.uRange]
w_range = [-1e-4, 8e-4]
# Pick cross section to view from file or default
manualMax = None

folders = args.files
labels = args.files
colors = ['xkcd:green','xkcd:red','xkcd:light blue','xkcd:blue']
lStyle = ['-','--']
resultFolder = '/results'
fileEnding = ""

# manualMax = 4957200

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

plotDPI = args.dpi
print('Plot DPI:',plotDPI)


figLabels = ["(a)","(b)","(c)","(d)","(e)","(f)"]
fig, axes = plt.subplots(2, 2, figsize=(10, 8), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[0,1]
# ax3 = axes[0,2]
ax4 = axes[1,0]
ax5 = axes[1,1]
# ax6 = axes[1,2]
#Open figure before looping
## Plotting U as function of depth
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
        if "dynDiag.0" in file:
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

    # maxFrames = args.numFrames[0]
    # if((maxStep-startStep)/sizeStep > maxFrames):   #if more than # frames, downscale to be less than #
    #     dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/maxFrames)
    #     print('Reducing time resolution by', dwnScale)
    #     sizeStep = sizeStep * dwnScale
    sizeStep = 57600 #force to 10 days time steps (assuming 15 dt in MITgcm)
    print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)
    
    if(manualMax != None):
        maxStep = manualMax #if manual adjusting needed
    z = np.squeeze(mds.rdmds("%s%s/RC" %(folder,resultFolder))) #midpoints of cells
    dz = np.load("%s/input/dz.npy" %(folder)) # we now save this upon grid generation, much better

    dynName = ['dynDiag']
    name = ['U']
    units = ["[m/s]"]

    stepsToAverage = 36
    tmpSize = np.array(np.shape(mds.rdmds("%s%s/dynDiag"%(folder,resultFolder), maxStep)))
    tmpSize = np.insert(tmpSize,0,0)
    tmpSize[0] = stepsToAverage
    print(tmpSize)
    dataArray = np.zeros(tmpSize)
    tmpSize[1] = 2
    dataArray_tracer = np.zeros(tmpSize)
    timeSteps = np.arange(maxStep-(stepsToAverage-1)*sizeStep, maxStep + 1, sizeStep)
    print("averaging over final",(timeSteps[-1]-timeSteps[0])*dt/86400,"days")
    for j in range(len(timeSteps)):
        dataArray[j,:,:,:,:] = mds.rdmds("%s%s/dynDiag"%(folder,resultFolder), timeSteps[j])[:,:,:,:]
        dataArray_tracer[j,:,:,:,:] = mds.rdmds("%s%s/ptraceDiag"%(folder,resultFolder), timeSteps[j])[:,:,:,:]
        dataArray[j,:,dataArray[j,1,:,:,:] == 0] = np.nan #cells in walls are set to NAN
        dataArray_tracer[j,:,dataArray[j,1,:,:,:] == 0] = np.nan #cells in walls are set to NAN
    MAM_data = np.nanmean(dataArray[0:9,:,:,:,:],axis=0)
    JJA_data = np.nanmean(dataArray[9:18,:,:,:,:],axis=0)
    SON_data = np.nanmean(dataArray[18:27,:,:,:,:],axis=0)
    DJF_data = np.nanmean(dataArray[27:36,:,:,:,:],axis=0)
    MAM_data_tracer = np.nanmean(dataArray_tracer[0:9,:,:,:,:],axis=0)
    JJA_data_tracer = np.nanmean(dataArray_tracer[9:18,:,:,:,:],axis=0)
    SON_data_tracer = np.nanmean(dataArray_tracer[18:27,:,:,:,:],axis=0)
    DJF_data_tracer = np.nanmean(dataArray_tracer[27:36,:,:,:,:],axis=0)
    # annual_data = np.nanmean(dataArray,axis=0)
    # mélange 1-10 km ()
    # ax1.plot(np.nanmean(MAM_data[0,:,:,1:25],axis=(1,2)),z,color=colors[0],linewidth=2,linestyle=lStyle[i])
    if(args.profile > 0  and i == 0):
        tt = np.load('2015_Melange_T.npy')
        zz = -1*np.load('2015_Melange_Z.npy')
        tt2 = np.load('2010mar_Melange_T.npy')
        zz2 = -1*np.load('2010mar_Z.npy')
        ax1.plot(tt,zz,color=colors[1],linewidth=2,linestyle=':')
        # ax1.plot(tt2,zz2,color=colors[3],linewidth=2,linestyle=':')
    ax1.plot(np.nanmean(JJA_data[0,:,:,1:25],axis=(1,2)),z,color=colors[1],linewidth=2,linestyle=lStyle[i])
    # ax1.plot(np.nanmean(SON_data[0,:,:,1:25],axis=(1,2)),z,color=colors[2],linewidth=2,linestyle=lStyle[i])
    ax1.plot(np.nanmean(DJF_data[0,:,:,1:25],axis=(1,2)),z,color=colors[3],linewidth=2,linestyle=lStyle[i])


    # ax2.plot(np.nanmean(MAM_data[2,:,:,1:25],axis=(1,2)),z,color=colors[0],linewidth=2,linestyle=lStyle[i])
    ax2.plot(np.nanmean(JJA_data[2,:,:,1:25],axis=(1,2)),z,color=colors[1],linewidth=2,linestyle=lStyle[i])
    # ax2.plot(np.nanmean(SON_data[2,:,:,1:25],axis=(1,2)),z,color=colors[2],linewidth=2,linestyle=lStyle[i])
    ax2.plot(np.nanmean(DJF_data[2,:,:,1:25],axis=(1,2)),z,color=colors[3],linewidth=2,linestyle=lStyle[i])

    # # ax3.plot(np.nanmean(MAM_data[3,:,:,1:25],axis=(1,2)),z,color=colors[0],linewidth=2,linestyle=lStyle[i])
    # ax3.plot(np.nanmean(JJA_data_tracer[0,:,:,1:25],axis=(1,2)),z,color=colors[1],linewidth=2,linestyle=lStyle[i])
    # # ax3.plot(np.nanmean(SON_data[3,:,:,1:25],axis=(1,2)),z,color=colors[2],linewidth=2,linestyle=lStyle[i])
    # ax3.plot(np.nanmean(DJF_data_tracer[0,:,:,1:25],axis=(1,2)),z,color=colors[3],linewidth=2,linestyle=lStyle[i])
    # # ax3.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    # mid fjord 30-60 km
    # ax4.plot(np.nanmean(MAM_data[0,:,:,75:150],axis=(1,2)),z,color=colors[0],linewidth=2,linestyle=lStyle[i])
    if(args.profile > 0 and i == 0):
        tt = np.load('2015_Fjord_T.npy')
        zz = -1*np.load('2015_Fjord_Z.npy')
        tt2 = np.load('2010mar_Fjord_T.npy')
        zz2 = -1*np.load('2010mar_Z.npy')
        ax4.plot(tt,zz,color=colors[1],linewidth=2,linestyle=':')
        # ax4.plot(tt2,zz2,color=colors[3],linewidth=2,linestyle=':')
    ax4.plot(np.nanmean(JJA_data[0,:,:,75:154],axis=(1,2)),z,color=colors[1],linewidth=2,linestyle=lStyle[i])
    # ax4.plot(np.nanmean(SON_data[0,:,:,75:150],axis=(1,2)),z,color=colors[2],linewidth=2,linestyle=lStyle[i])
    ax4.plot(np.nanmean(DJF_data[0,:,:,75:150],axis=(1,2)),z,color=colors[3],linewidth=2,linestyle=lStyle[i])

    if(args.labels != None):
        if(i==0):
            ax4.plot([],[],color='gray',label='xCTD',linestyle=':')
        ax4.plot([],[],color='gray',label=args.labels[i],linestyle=lStyle[i])

    # ax5.plot(np.nanmean(MAM_data[2,:,:,75:150],axis=(1,2)),z,color=colors[0],linewidth=2,linestyle=lStyle[i])
    ax5.plot(np.nanmean(JJA_data[2,:,:,75:150],axis=(1,2)),z,color=colors[1],linewidth=2,linestyle=lStyle[i])
    # ax5.plot(np.nanmean(SON_data[2,:,:,75:150],axis=(1,2)),z,color=colors[2],linewidth=2,linestyle=lStyle[i])
    ax5.plot(np.nanmean(DJF_data[2,:,:,75:150],axis=(1,2)),z,color=colors[3],linewidth=2,linestyle=lStyle[i])

    # # ax6.plot(np.nanmean(MAM_data_tracer[0,:,:,75:150],axis=(1,2)),z,color=colors[0],linewidth=2,linestyle=lStyle[i])
    # ax6.plot(np.nanmean(JJA_data_tracer[0,:,:,75:150],axis=(1,2)),z,color=colors[1],linewidth=2,linestyle=lStyle[i])
    # # ax6.plot(np.nanmean(SON_data_tracer[0,:,:,75:150],axis=(1,2)),z,color=colors[2],linewidth=2,linestyle=lStyle[i])
    # ax6.plot(np.nanmean(DJF_data_tracer[0,:,:,75:150],axis=(1,2)),z,color=colors[3],linewidth=2,linestyle=lStyle[i])

# ax1.plot([],[],label='MAM',linestyle='-',color=colors[0])
ax1.plot([],[],label='Summer (JJA)',linestyle='-',color=colors[1])
# ax1.plot([],[],label='SON',linestyle='-',color=colors[2])
ax1.plot([],[],label='Winter (DJF)',linestyle='-',color=colors[3])

yRange = [-405,10]

ax1.grid(alpha=.5)
top_title = 'Mélange [0-10 km]'
bottom_title = 'Mid Fjord [30-60 km]'
ax1.set_title(top_title)
ax1.set_ylabel('Depth [m]')
ax1.set_xlabel('Temp [C]')
ax1.set_xlim(t_range)
ax1.set_ylim(yRange)
ax1.legend()
ax2.grid(alpha=.5)
ax2.set_title(top_title)
ax2.set_ylabel('Depth [m]')
ax2.set_xlabel('U [m/s]')
ax2.set_xlim(u_range)
ax2.set_ylim(yRange)
# ax3.grid(alpha=.5)
# ax3.set_title(top_title)
# ax3.set_ylabel('Depth [m]')
# ax3.set_xlabel('SGD Fraction [ ]')
# ax3.set_xlim([-.01, .04])
ax4.grid(alpha=.5)
ax4.set_title(bottom_title)
ax4.set_ylabel('Depth [m]')
ax4.set_xlabel('Temp [C]')
ax4.set_xlim(t_range)
ax4.set_ylim(yRange)
if(args.labels != None):
    ax4.legend()
ax5.grid(alpha=.5)
ax5.set_title(bottom_title)
ax5.set_ylabel('Depth [m]')
ax5.set_xlabel('U [m/s]')
ax5.set_xlim(u_range)
ax5.set_ylim(yRange)
# ax6.grid(alpha=.5)
# ax6.set_title(bottom_title)
# ax6.set_ylabel('Depth [m]')
# ax6.set_xlabel('SGD Fraction [ ]')
# ax6.set_xlim([-.01, .04])

# for label,ax in zip(figLabels,axes.flatten()):
#     ax.text(
#         ax.get_xlim()[0], ax.get_ylim()[1], label,
#         fontsize='x-large', va='bottom',ha='center', fontfamily='sans serif')
    
timeString = ''
if(args.timeRange != None):
    timeString = f'_T_{args.timeRange[0]}_{args.timeRange[1]}'

plt.savefig(f'figs/FigureOSM{fileEnding}{timeString}.png', format='png',dpi=plotDPI)
if(args.silent == 0):
    plt.show()
plt.close()
   