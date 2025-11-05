#!/usr/bin/env python3

# Paul Summers November 2025
# Made for making figures for melange collapse paper


import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import sys
from scipy.integrate import simpson
import argparse
import os
import fileinput
# sys.path.append('/hdd/glaciome/models/glaciome1D')
# sys.path.insert(0, '')
sys.path.append('/Users/psummers8/Documents/glaciome1D')
from glaciome1D import constants, glaciome

import glob
import pickle

saveProfiles = False

constant = constants()

# Take input options. Some defaults are set here, so be aware
parser = argparse.ArgumentParser(description='Plot options for melange changed over time')
parser.add_argument('-f','--files', nargs='*', default=['couplingResults/MITgcmRun_'],
                    help='file string used to name pickle files FILES00000.pickle [default = couplingResults/MITgcmRun_]')
parser.add_argument('-l','--labels', nargs='*', default=None,
                    help='file string label plots [default = None]')
parser.add_argument('-sgd','--sgd', nargs='*', type=float, default = None,
                    help='SGD, any value loads SGD from input, needs results RY')
parser.add_argument('-s','--silent', action='count', default=0,
                    help='Option to silence showing of plots')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = [0,  0],
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('-w','--windowMean', nargs='?', type=int, default=10,
                    help='window width for time averaging melt rates [default = 10]')
parser.add_argument('-ff','--fixForce', action='count', default=0,
                    help='Option to fix forcing left scale to +/- 5%%')
parser.add_argument('-ll','--legLoc', nargs='?', default=0,type=int,
                    help='legend location [default = 0]')
parser.add_argument('-dpi','--dpi', nargs='?', default=200,type=int,
                    help='dpi to print [default = 200]')
args = parser.parse_args()

# print(args)
figLabels = ["(a)","(b)","(c)","(d)","(e)","(f)"]
fig, axes = plt.subplots(2, 2, figsize=(10, 6), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[1,0]
ax3 = axes[0,1]
ax4 = axes[1,1]

lStyle = ['-','--',':','-.',(0, (3, 2, 1, 2, 1, 2)),(0, (2, 3, 1, 2, 1, 2, 1, 2,))]
meltColors = ['xkcd:tomato red','xkcd:rose','xkcd:wine','xkcd:grape purple','xkcd:grape','xkcd:violet']

colorList = ['xkcd:light blue',
                'xkcd:bright orange','xkcd:tomato',
                'xkcd:faded pink','xkcd:dark blue']

if(args.silent > 0): #print status for script
    print(f"rampPlot: {args}")

for i in range(len(args.files)):
    fileStr = args.files[i]
    files = sorted(glob.glob('%s*.pickle'%fileStr))
    jShift = 0
    if(args.timeRange == None):
        toIterate = np.arange(len(files))
    else:
        if(args.timeRange[1] == 0):
            toIterate = np.arange(len(files[args.timeRange[0]:]))
        else:
            toIterate = np.arange(len(files[args.timeRange[0]:args.timeRange[1]]))
        jShift = args.timeRange[0]
    H0Time = np.zeros(np.shape(toIterate))
    VTime = np.zeros(np.shape(toIterate))
    UcTime = np.zeros(np.shape(toIterate))
    UtTime = np.zeros(np.shape(toIterate))
    lengthTime = np.zeros(np.shape(toIterate))
    iterationNumber = np.zeros(np.shape(toIterate))
    pressTime = np.zeros(np.shape(toIterate))
    BTime = np.zeros(np.shape(toIterate))
    BShapeTime = np.zeros([np.shape(toIterate)[0],20])
    gTime = np.zeros(np.shape(toIterate))
    spdTime = np.zeros(np.shape(toIterate))
    timeTime = np.zeros(np.shape(toIterate))
    muSTime = np.zeros(np.shape(toIterate))
    xTermTime = np.zeros(np.shape(toIterate))
    bFluxTime = np.zeros(np.shape(toIterate))
    endFluxTime = np.zeros(np.shape(toIterate))
    inFluxTime = np.zeros(np.shape(toIterate))
    for j in toIterate:
        file = files[j + jShift]
        with open(files[j + jShift], 'rb') as fileName:
            data = pickle.load(fileName)
            fileName.close()
        if(j == toIterate[-1]):
            print(f"Final state: {data}")
        it = file.replace(fileStr,'').replace('./', '').replace('.pickle', '')
        H0Time[j]=data.H0
        xTermTime[j] = data.X[0]
        X_ = np.concatenate(([data.X[0]], data.X_, [data.X[-1]]))
        H = np.concatenate(([data.H0], data.H, [data.HL]))
        W = np.concatenate(([data.W0], data.W, [data.WL]))
        VTime[j] = simpson(H*W, x=X_)*1e-9
        lengthTime[j]=data.X[-1] - data.X[0]  #difference between these
        UcTime[j] = data.Uc
        UtTime[j] = data.Ut
        BTime[j] = -1 * np.mean(data.B) /365.25 # we flip this for plotting purposes
        BShapeTime[j,:] = -1 * data.B /365.25
        timeTime[j] = data.t*365 
        pressTime[j] = data.force()#data.H0*data.pressure(data.H0) #should I use force?
        spdTime[j]=np.mean(data.U)
        gTime[j]=np.mean(data.gg)
        muSTime[j] = data.param.muS
        iterationNumber[j]=int(it)
        endFluxTime[j] = data.HL*data.WL*data.U[-1] / 3.154e7 #m^3/s
        if(j != toIterate[0]): #include mass loss to length loss
            endFluxTime[j] += (data.H[-2])*data.WL * (lengthTime[j-1] - lengthTime[j]) / (data.dt * 3.154e7) #m^3/s
        inFluxTime[j] = data.Ht*data.Uc*data.WL / 3.154e7 #m^3/s
        bFluxTime[j] = -1*simpson(data.B*data.W, x=data.X_)/ 3.154e7  #m^3/s
    timeTime = timeTime - timeTime[0]
    # if(not args.timeRange == None):
    #     timeTime = timeTime + args.timeRange[0]
    timeLabel = 'Time [days]'
    
    from scipy import interpolate
    from MITgcmutils import mds
    remoteResult = args.files[i].rsplit('/',1)[0]
    print(f'remoteResult is ({remoteResult})')
    y = np.squeeze(mds.rdmds(f"{remoteResult}/../results/YC")[:,0])
    dy = 2*y[0]
    for line in fileinput.input(f'{remoteResult}/../input/data'):
        if "ExternForcingPeriod=" in line:
            period = float(line[21:-2])
        elif "ExternForcingCycle=" in line:
            cycle = float(line[20:-2])
    nt = int(cycle/period)
    plumeMask = np.fromfile(f'{remoteResult}/../input/plumeMask.bin', dtype='>f8') #(2 is line, 3 is semi-cone)
    rad = np.fromfile(f'{remoteResult}/../input/runoffRad.bin', dtype='>f8')
    rad = rad.reshape((nt,len(plumeMask)))
    vel = np.fromfile(f'{remoteResult}/../input/runoffVel.bin', dtype='>f8')
    vel = vel.reshape((nt,len(plumeMask)))

    plumeType = np.max(plumeMask)
    cleanRad = rad[:,plumeMask == plumeType]
    cleanVel = vel[:,plumeMask == plumeType]
    seasonTime = np.linspace(0,cycle/86400,nt)
    if(plumeType == 2):
        runoff = dy*cleanRad[:,0]*cleanVel[:,0]
    else: #semi-rad
        runoff = np.pi*cleanRad**2*cleanVel
    # print(f'seasonTime {seasonTime}, nt {nt}')
    # print(f'runoff {runoff}')
    # print(f'plume type {plumeType}')
    # print(np.shape(seasonTime))
    # print(np.shape(runoff))
    f = interpolate.interp1d(seasonTime, runoff,fill_value='0')
    if(args.sgd != None):
        timeTime = f(timeTime % (cycle/86400))
        timeLabel = 'SGD [$m^3/s$]'

    N = args.windowMean #window over to smooth in time points (days)
    #Pad edges of BTime to have smooth edges and work for short time series
    padL = int(N/2)
    padR = int((N-1)/2)  #This ensures odds work OK
    BTimeSmooth = np.convolve(np.concatenate((BTime[0]*np.ones(padL),BTime,BTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    # timeTime = (timeTime - timeTime[0])
    forceTimeRough = BTime
    forceTime = BTimeSmooth
    bFluxSmooth = np.convolve(np.concatenate((bFluxTime[0]*np.ones(padL),bFluxTime,bFluxTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    endFluxSmooth = np.convolve(np.concatenate((endFluxTime[0]*np.ones(padL),endFluxTime,endFluxTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')

    forceLabel = 'Melt Rate [m/day]'
    forceString = 'Melt'
    notforceTime = UcTime
    notforceLabel = 'Calving Rate [m/yr]'
    notforceString = 'Calving'

    # plt.suptitle('Mélange over Time')
    if(args.labels == None):
        ax1.plot(timeTime,forceTime,color=meltColors[i],linestyle=lStyle[i])
    else: #labels names have been supplied
        ax1.plot(timeTime,forceTime,color=meltColors[0],linestyle=lStyle[i],label=args.labels[i])
    ax1.plot(timeTime,forceTimeRough,color=meltColors[0],linestyle='-',alpha=.25)
    ax1.scatter(timeTime[0],forceTime[0],s=50,marker='*',color='black')
    
    ax1.set_ylabel(forceLabel,color=meltColors[0])
    ax1.tick_params(axis='y',labelcolor=meltColors[0])
    if(args.fixForce > 0):
        ax1.set_ylim([0.95 * np.nanmean(forceTimeRough), 1.05* np.nanmean(forceTimeRough)])
    ax1.set_xlabel(timeLabel)
    ax1.grid(alpha=.5)
    ax1.set_xticks([0,73,73*2,73*3,73*4,73*5])
    ax1_2 = ax1.twinx()
    ax1_2.plot(timeTime,f((timeTime  + args.timeRange[0]) % (cycle/86400)),color='xkcd:azure',linestyle=lStyle[i],alpha=.5)
    ax1_2.set_ylabel('SGD [$m^3/s$]',color='xkcd:azure')
    ax1_2.tick_params(axis='y',labelcolor='xkcd:azure')

    # ax1.set_title('Melt')
    if(args.labels != None):
        ax1.legend(loc=args.legLoc)
    xString =''

    drivingVariable = timeTime
    drivingLabel = timeLabel
    xString = 'Time'
    if(args.sgd != None):
        xString = 'SGD'
    roughLT = lengthTime
    roughHT = H0Time
    lengthTime = np.convolve(np.concatenate((lengthTime[0]*np.ones(padL),lengthTime,lengthTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    H0Time = np.convolve(np.concatenate((H0Time[0]*np.ones(padL),H0Time,H0Time[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    ax2.plot(drivingVariable,roughLT,color='xkcd:blueberry',alpha=0.25)
    ax2.plot(drivingVariable,lengthTime,color='xkcd:blueberry',linestyle=lStyle[i])
    ax2.scatter(drivingVariable[0],lengthTime[0],s=50,marker='*',color='black')
    ax2.set_ylabel('Mélange length [m]',color='xkcd:blueberry')
    ax2.tick_params(axis='y',labelcolor='xkcd:blueberry')
    ax2.grid(alpha=.5)
    ax2.set_xticks([0,73,73*2,73*3,73*4,73*5])
    # sca=ax2.scatter(BTime,lengthTime,s=None,c=timeTime,cmap='cividis')
    # cbar=plt.colorbar(sca)
    # cbar.set_label('Iteration')
    if(True): #gate patch addition
        from matplotlib import patches
        init_xlim = ax2.get_xlim()
        init_ylim = ax2.get_ylim()
        p1 = patches.Rectangle([0   ,init_ylim[0]],73,init_ylim[1]-init_ylim[0],facecolor=colorList[0],alpha=.1,edgecolor=None)
        p2 = patches.Rectangle([73  ,init_ylim[0]],73,init_ylim[1]-init_ylim[0],facecolor=colorList[1],alpha=.1,edgecolor=None)
        p3 = patches.Rectangle([73*2,init_ylim[0]],73,init_ylim[1]-init_ylim[0],facecolor=colorList[2],alpha=.1,edgecolor=None)
        p4 = patches.Rectangle([73*3,init_ylim[0]],73,init_ylim[1]-init_ylim[0],facecolor=colorList[3],alpha=.1,edgecolor=None)
        p5 = patches.Rectangle([73*4,init_ylim[0]],73,init_ylim[1]-init_ylim[0],facecolor=colorList[4],alpha=.1,edgecolor=None)
        ax2.add_patch(p1)
        ax2.add_patch(p2)
        ax2.add_patch(p3)
        ax2.add_patch(p4)
        ax2.add_patch(p5)




    if(False): #gate for plotting thickness
        if(i == 0):
            ax2_2 = ax2.twinx()
        ax2_2.plot(drivingVariable,roughHT,color='xkcd:mulberry',alpha=0.25)
        ax2_2.plot(drivingVariable,H0Time,color='xkcd:mulberry',linestyle=lStyle[i])
        ax2_2.scatter(drivingVariable[0],H0Time[0],s=50,marker='*',color='black')
        ax2_2.set_ylabel('Mélange H0 [m]',color='xkcd:mulberry')
        ax2_2.tick_params(axis='y',labelcolor='xkcd:mulberry')
        # ax2_2.grid(alpha=.25,color='xkcd:mulberry')
    ax2.set_xlabel(drivingLabel)
    # ax2.set_title('Size')

    if(i == 0): #Plot blanks for legend 
        ax3.plot([],[],color='xkcd:sand',linestyle='-',label='In')   
        ax3.plot([],[],color='xkcd:apple',linestyle='-',label='Out')  
        ax3.plot([],[],color='xkcd:red',linestyle='-',label='Melt')  
        ax3.plot([],[],color='xkcd:sky blue',linestyle='-',label='Net',alpha=.5)
        # ax3.plot([],[],color='xkcd:indigo',linestyle='-',label='Volume')
    ax3.plot(timeTime,inFluxTime,color='xkcd:sand',linestyle=lStyle[i])   
    ax3.plot(timeTime,endFluxTime,color='xkcd:apple',linestyle=lStyle[i],alpha = .25)  
    ax3.plot(timeTime,bFluxTime,color='xkcd:red',linestyle=lStyle[i],alpha = .25)   
    ax3.plot(timeTime,endFluxSmooth,color='xkcd:apple',linestyle=lStyle[i])  
    ax3.plot(timeTime,bFluxSmooth,color='xkcd:red',linestyle=lStyle[i])
    ax3.plot(timeTime,inFluxTime-endFluxSmooth-bFluxSmooth,color='xkcd:sky blue',linestyle=lStyle[i],alpha=.75)
    ax3.plot(timeTime,inFluxTime-endFluxTime-bFluxTime,color='xkcd:sky blue',linestyle=lStyle[i],alpha=.25)
    ax3.plot([timeTime[0],timeTime[-1]],[0,0],color = 'gray',alpha=.5,linestyle='--')    
    ax3.set_ylabel('Volume Flux [m^3/s]')
    ax3.set_xlabel('Time [days]')
    ax3.grid(alpha=.5)
    # ax3.set_title('')
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=4)
    ax3.set_xticks([0,73,73*2,73*3,73*4,73*5])
    
    # for ii in range(len(timeTime)):
    #     ax4.plot(np.linspace(0,1,20),BShapeTime[ii,:],alpha = .1,color=colorList[int(timeTime[ii]/73)%5])
    ax4.plot(np.linspace(0,1,20),np.mean(BShapeTime[  0: 73,:],axis=0),alpha = 1,color=colorList[0],label='Late Winter')
    ax4.plot(np.linspace(0,1,20),np.mean(BShapeTime[ 73:164,:],axis=0),alpha = 1,color=colorList[1],label='Plume Rise')
    ax4.plot(np.linspace(0,1,20),np.mean(BShapeTime[164:219,:],axis=0),alpha = 1,color=colorList[2],label='Plume Peak')
    ax4.plot(np.linspace(0,1,20),np.mean(BShapeTime[219:292,:],axis=0),alpha = 1,color=colorList[3],label='Plume Fall')
    ax4.plot(np.linspace(0,1,20),np.mean(BShapeTime[292:365,:],axis=0),alpha = 1,color=colorList[4],label='Early Winter')
    ax4.grid(alpha=.5)
    ax4.legend(loc=9)
    # if(i == 0):
    #     ax4_2=ax4.twinx()
    # ax4_2.plot(drivingVariable,roughG,color='xkcd:lilac',alpha=0.25)
    # ax4_2.plot(drivingVariable,gTime,color='xkcd:lilac',linestyle=lStyle[i])
    # ax4_2.scatter(drivingVariable[0],gTime[0],s=50,marker='*',color='black')
    # ax4_2.set_ylabel('fluidity $g^{\\prime}$ [$yr^{-1}$]   ',color='xkcd:lilac')
    # ax4_2.tick_params(axis='y',labelcolor='xkcd:lilac')
    # ax4.set_title('Avg Speed and Fluidity $g^{\\prime}$')
    ax4.set_xlabel('Fraction along mélange [ ]')
    ax4.set_ylabel('Melt Rate [m/day]')
    ax4.set_ylim([0,1.6])

for label,ax in zip(figLabels,axes.reshape(4,1)):
    ax=ax[0]
    ax.text(
        ax.get_xlim()[0], ax.get_ylim()[1], label,
        fontsize='x-large', va='bottom',ha='right', fontfamily='sans serif')

dirStr = ''
if(os.path.isdir('figs')):
    dirStr = 'figs/'
multiStr = ''
timeString = ''
if(i > 0):
    multiStr = f'Multi{i+1}'
if(args.timeRange != None):
    multiStr = f'{multiStr}_T_{args.timeRange[0]}_{args.timeRange[1]}'

plt.savefig('%sprettyPlot%s%s%s.png' %(dirStr,forceString,xString,multiStr),format='png',dpi=args.dpi)
if(args.silent == 0):
    plt.show()
plt.close()
