#!/usr/bin/env python3

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
parser.add_argument('-x','--xVariable', nargs=1, type=int, default=[1],
                    help='what to plot against 0: xVar , 1: Time [default]')
parser.add_argument('-sgd','--sgd', nargs='*', type=float, default = None,
                    help='SGD, any value loads SGD from input, needs results RY')
parser.add_argument('-r','--rampVariable', nargs=1, type=int, default=[0],
                    help='what is ramping 0: MeltRate [default], 1: Uc')
parser.add_argument('-s','--silent', action='count', default=0,
                    help='Option to silence showing of plots')
parser.add_argument('-ms','--muS', action='count', default=0,
                    help='Option to show muS on forcing plots')
parser.add_argument('-xt','--xTerminus', action='count', default=0,
                    help='Option to show terminus location on plots')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
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
meltColors = ['xkcd:tomato','xkcd:rose','xkcd:wine','xkcd:grape purple','xkcd:grape','xkcd:violet']

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
    gTime = np.zeros(np.shape(toIterate))
    spdTime = np.zeros(np.shape(toIterate))
    timeTime = np.zeros(np.shape(toIterate))
    muSTime = np.zeros(np.shape(toIterate))
    xTermTime = np.zeros(np.shape(toIterate))
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
        BTime[j] = -1 * np.mean(data.B) /365.0 # we flip this for plotting purposes
        timeTime[j] = data.t*365.0 
        pressTime[j] = data.force()#data.H0*data.pressure(data.H0) #should I use force?
        spdTime[j]=np.mean(data.U)
        gTime[j]=np.mean(data.gg)
        muSTime[j] = data.param.muS
        iterationNumber[j]=int(it)
    timeTime = timeTime - timeTime[0] + jShift# shift to start at t=0

    timeLabel = 'Time [days]'
    if(args.sgd != None):
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
        # print(f'seasonTime {seasonTime}, nt {nt}')
        # print(f'runoff {runoff}')
        # print(f'plume type {plumeType}')
        # print(np.shape(seasonTime))
        # print(np.shape(runoff))
        f = interpolate.interp1d(seasonTime, runoff,fill_value='0')
        timeTime = f(timeTime)
        timeLabel = 'SGD [$m^3/s$]'

    N = args.windowMean #window over to smooth in time points (days)
    #Pad edges of BTime to have smooth edges and work for short time series
    padL = int(N/2)
    padR = int((N-1)/2)  #This ensures odds work OK
    BTimeSmooth = np.convolve(np.concatenate((BTime[0]*np.ones(padL),BTime,BTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    # timeTime = (timeTime - timeTime[0])/5 #if pre-code fix for total time
    if(args.rampVariable[0] == 0):
        forceTimeRough = BTime
        forceTime = BTimeSmooth
        forceLabel = 'Melt Rate [m/day]'
        forceString = 'Melt'
        notforceTime = UcTime
        notforceLabel = 'Calving Rate [m/yr]'
        notforceString = 'Calving'
    elif(args.rampVariable[0] == 1):
        forceTime = UcTime
        forceLabel = 'Calving Rate [m/yr]'
        forceString = 'Calving'
        notforceTime = BTimeSmooth
        notforceLabel = 'Melt Rate [m/day]'
        notforceString = 'Melt'
    else:
        raise Exception('Invalid forcing option %i' %args.deltaVariable[0])

    plt.suptitle('Mélange with Time Varying %s' %forceString)
    if(args.labels == None):
        ax1.plot(timeTime,forceTime,color=meltColors[i],linestyle=lStyle[i])
    else: #labels names have been supplied
        ax1.plot(timeTime,forceTime,color=meltColors[i],linestyle=lStyle[i],label=args.labels[i])
    if(args.rampVariable[0] == 0):
        ax1.plot(timeTime,forceTimeRough,color=meltColors[i],linestyle='-',alpha=.25)
    ax1.scatter(timeTime[0],forceTime[0],s=50,marker='*',color='black')
    # sca=ax1.scatter(BTime,H0Time,s=None,c=timeTime,cmap='cividis')
    # cbar=plt.colorbar(sca)
    # cbar.set_label('Iteration')
    if(args.muS > 0):
        if(i == 0):
            ax1_2 = ax1.twinx()
        ax1_2.plot(timeTime,muSTime,color='xkcd:azure',linestyle=lStyle[i])
        ax1_2.scatter(timeTime[0],muSTime[0],s=50,marker='*',color='black')
        ax1_2.set_ylabel('$\\mu_S$ [ ]',color='xkcd:azure')
        ax1_2.tick_params(axis='y',labelcolor='xkcd:azure')
        # ax1_2.grid(alpha=.25,color='xkcd:azure')
        
    #If second axis is free and non-forced variable (B or Uc) changes, plot it
    elif(True):
        if(i == 0):
            ax1_2 = ax1.twinx()
        roughNFT = notforceTime
        notforceTime = np.convolve(np.concatenate((notforceTime[0]*np.ones(padL),notforceTime,notforceTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
        ax1_2.plot(timeTime,roughNFT,color='xkcd:azure',alpha=0.25)
        ax1_2.plot(timeTime,notforceTime,color='xkcd:azure',linestyle=lStyle[i])
        ax1_2.scatter(timeTime[0],notforceTime[0],s=50,marker='*',color='black')
        ax1_2.set_ylabel(notforceLabel,color='xkcd:azure')
        ax1_2.tick_params(axis='y',labelcolor='xkcd:azure')
        if(args.xTerminus > 0 and i == 0):
            ax1_2.plot(timeTime,UtTime,color='xkcd:gray',linestyle=lStyle[i])
            ax1_2.scatter(timeTime[0],UtTime[0],s=50,marker='*',color='black')
        # ax1_2.set_ylim([5900, 8800])
    ax1.set_ylabel(forceLabel)
    if(args.fixForce > 0):
        ax1.set_ylim([0.95 * np.nanmean(forceTimeRough), 1.05* np.nanmean(forceTimeRough)])
    ax1.set_xlabel(timeLabel)
    ax1.grid(alpha=.5)
    ax1.set_title('%s Forcing' %forceString)
    if(args.labels != None):
        ax1.legend(loc=args.legLoc)
    xString =''
    if(args.xVariable[0] == 0):
        drivingVariable = forceTime
        drivingLabel = forceLabel
        xString = forceString
    elif(args.xVariable[0] == 1):
        drivingVariable = timeTime
        drivingLabel = timeLabel
        xString = 'Time'
        if(args.sgd != None):
            xString = 'SGD'
    else:
        raise Exception('Invalid X-axis option: %s' %args.xVariable[0])
    roughLT = lengthTime
    roughHT = H0Time
    lengthTime = np.convolve(np.concatenate((lengthTime[0]*np.ones(padL),lengthTime,lengthTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    H0Time = np.convolve(np.concatenate((H0Time[0]*np.ones(padL),H0Time,H0Time[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    ax2.plot(drivingVariable,roughLT,color='xkcd:blueberry',alpha=0.25)
    ax2.plot(drivingVariable,lengthTime,color='xkcd:blueberry',linestyle=lStyle[i])
    ax2.scatter(drivingVariable[0],lengthTime[0],s=50,marker='*',color='black')
    ax2.set_ylabel('Mélange length [m]',color='xkcd:blueberry')
    ax2.tick_params(axis='y',labelcolor='xkcd:blueberry')
    ax2.grid(alpha=.25,color='xkcd:blueberry')
    # sca=ax2.scatter(BTime,lengthTime,s=None,c=timeTime,cmap='cividis')
    # cbar=plt.colorbar(sca)
    # cbar.set_label('Iteration')
    if(i == 0):
        ax2_2 = ax2.twinx()
    ax2_2.plot(drivingVariable,roughHT,color='xkcd:mulberry',alpha=0.25)
    ax2_2.plot(drivingVariable,H0Time,color='xkcd:mulberry',linestyle=lStyle[i])
    ax2_2.scatter(drivingVariable[0],H0Time[0],s=50,marker='*',color='black')
    ax2_2.set_ylabel('Mélange H0 [m]',color='xkcd:mulberry')
    ax2_2.tick_params(axis='y',labelcolor='xkcd:mulberry')
    # ax2_2.grid(alpha=.25,color='xkcd:mulberry')
    ax2.set_xlabel(drivingLabel)
    ax2.set_title('Size')

    roughPress = pressTime
    pressTime = np.convolve(np.concatenate((pressTime[0]*np.ones(padL),pressTime,pressTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    ax3.plot(drivingVariable,roughPress,color='xkcd:sunflower',alpha=0.25)
    ax3.plot(drivingVariable,pressTime,color='xkcd:sunflower',linestyle=lStyle[i])
    ax3.scatter(drivingVariable[0],pressTime[0],s=50,marker='*',color='black')
    # sca=ax3.scatter(BTime,pressTime,s=None,c=timeTime,cmap='cividis')
    # cbar=plt.colorbar(sca)
    # cbar.set_label('Time [days]')
    ax3.set_ylabel('Mélange F/W [$Nm^{-1}$]')
    # ax3.set_yscale('log')
    ax3.set_xlabel(drivingLabel)
    ax3.grid(alpha=.5)
    ax3.set_title('Buttressing Strength')
    if(args.xTerminus > 0):
        if(i == 0):
            ax3_2 = ax3.twinx()
        ax3_2.plot(timeTime,xTermTime,color='xkcd:pumpkin',linestyle=lStyle[i])
        ax3_2.scatter(timeTime[0],xTermTime[0],s=50,marker='*',color='black')
        ax3_2.set_ylabel('Terminus position [m]',color='xkcd:pumpkin')
        ax3_2.tick_params(axis='y',labelcolor='xkcd:pumpkin')
        # ax1_2.grid(alpha=.25,color='xkcd:azure')

    roughSPD = spdTime
    roughG = gTime
    spdTime = np.convolve(np.concatenate((spdTime[0]*np.ones(padL),spdTime,spdTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    gTime = np.convolve(np.concatenate((gTime[0]*np.ones(padL),gTime,gTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    ax4.plot(drivingVariable,roughSPD/365.0,color='xkcd:apple',alpha=0.25)
    ax4.plot(drivingVariable,spdTime/365.0,color='xkcd:apple',linestyle=lStyle[i])
    ax4.scatter(drivingVariable[0],spdTime[0]/365.0,s=50,marker='*',color='black')
    ax4.set_ylabel('Speed [m/day]',color='xkcd:apple')
    ax4.tick_params(axis='y',labelcolor='xkcd:apple')
    ax4.grid(alpha=.25,color='xkcd:apple')
    if(i == 0):
        ax4_2=ax4.twinx()
    # sca=ax4.scatter(spdTime/365.0,gTime,s=20,c=timeTime,cmap='cividis')
    # cbar=plt.colorbar(sca)
    # cbar.set_label('Time [days]')
    ax4_2.plot(drivingVariable,roughG,color='xkcd:lilac',alpha=0.25)
    ax4_2.plot(drivingVariable,gTime,color='xkcd:lilac',linestyle=lStyle[i])
    ax4_2.scatter(drivingVariable[0],gTime[0],s=50,marker='*',color='black')
    ax4_2.set_ylabel('$g^{\\prime}$    ',color='xkcd:lilac')
    ax4_2.tick_params(axis='y',labelcolor='xkcd:lilac')
    # ax4_2.grid(alpha=.25,color='xkcd:lilac')
    ax4.set_title('Avg Speed and Fluidity $g^{\\prime}$')
    ax4.set_xlabel(drivingLabel)

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

plt.savefig('%sautoRamp%sResults%s%s.png' %(dirStr,forceString,xString,multiStr),format='png',dpi=args.dpi)
if(args.silent == 0):
    plt.show()
plt.close()
