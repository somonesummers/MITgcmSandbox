#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import sys
from scipy.integrate import simpson
import argparse
import os
# sys.path.append('/hdd/glaciome/models/glaciome1D')
# sys.path.insert(0, '')
sys.path.append('/Users/psummers8/Documents/glaciome1D')
from glaciome1D import constants, glaciome

import glob
import pickle

saveProfiles = False

constant = constants()

# Take input options. Some defaults are set here, so be aware
parser = argparse.ArgumentParser(description='Plot options for melange mass changed over time')
parser.add_argument('-f','--files', nargs='*', default=['couplingResults/MITgcmRun_'],
                    help='file string used to name pickle files FILES00000.pickle [default = couplingResults/MITgcmRun_]')
parser.add_argument('-l','--labels', nargs='*', default=None,
                    help='file string label plots [default = None]')
parser.add_argument('-s','--silent', action='count', default=0,
                    help='Option to silence showing of plots')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('-w','--windowMean', nargs='*', type=int, default=10,
                    help='window width for time averaging melt rates [default = 10]')
args = parser.parse_args()

# print(args)

fig, axes = plt.subplots(2, 2, figsize=(10, 6), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[1,0]
ax3 = axes[0,1]
ax4 = axes[1,1]

lStyle = ['-','--',':','-.']

if(args.silent > 0): #print status for script
    print(f"rampPlot: {args}")
showVolume = True

for i in range(len(args.files)):
    fileStr = args.files[i]
    files = sorted(glob.glob('%s*.pickle'%fileStr))
    jShift = 0
    if(args.timeRange == None):
        toIterate = np.arange(len(files))
    else:
        toIterate = np.arange(len(files[args.timeRange[0]:args.timeRange[1]]))
        jShift = args.timeRange[0]
    VTime = np.zeros(np.shape(toIterate))
    iterationNumber = np.zeros(np.shape(toIterate))
    bFluxTime = np.zeros(np.shape(toIterate))
    timeTime = np.zeros(np.shape(toIterate))
    endFluxTime = np.zeros(np.shape(toIterate))
    inFluxTime = np.zeros(np.shape(toIterate))
    H0Time = np.zeros(np.shape(toIterate))
    UcTime = np.zeros(np.shape(toIterate))
    lengthTime = np.zeros(np.shape(toIterate))
    BTime = np.zeros(np.shape(toIterate))
    pressTime = np.zeros(np.shape(toIterate))
    for j in toIterate:
        file = files[j + jShift]
        with open(files[j + jShift], 'rb') as fileName:
            data = pickle.load(fileName)
            fileName.close()
        if(j == toIterate[-1]):
            print(f"Final state: {data}")
        it = file.replace(fileStr,'').replace('./', '').replace('.pickle', '')
        H = np.concatenate(([data.H0], data.H, [data.HL]))
        W = np.concatenate(([data.W0], data.W, [data.WL]))
        VTime[j] = simpson(data.H*data.W, x=data.X_) #m^3
        endFluxTime[j] = data.HL*data.WL*data.U[-1] / 3.154e7 #m^3/s
        inFluxTime[j] = data.Ht*data.Uc*data.WL / 3.154e7 #m^3/s
        bFluxTime[j] = -1*simpson(data.B*data.W, x=data.X_)/ 3.154e7  #m^3/s
        BTime[j] = -1 * np.mean(data.B) /365.0 # we flip this for plotting purposes
        H0Time[j]=data.H0
        lengthTime[j]=data.X[-1] - data.X[0]  #difference between these
        UcTime[j] = data.Uc
        timeTime[j] = data.t*365.0 
        pressTime[j] = data.force()#data.H0*data.pressure(data.H0) #should I use force?

        iterationNumber[j]=int(it)
    timeTime = timeTime - timeTime[0] + jShift# shift to start at t=0# shift to start at t=0

    N = args.windowMean #window over to smooth in time points (days)
    BTimeSmooth = np.convolve(BTime, np.ones(N)/N, mode='same')
    bFluxTime = np.convolve(bFluxTime, np.ones(N)/N, mode='same')
    plt.suptitle('Mélange Volume Over Time')

    if(args.labels != None):
        lbTemp = args.labels[i]
    else:
        lbTemp = ''
    ax1.plot(timeTime,BTime,color='xkcd:red',linestyle='-',alpha = .25)
    ax1.plot(timeTime,BTimeSmooth,color='xkcd:red',linestyle=lStyle[i],label=lbTemp)
    # ax1.scatter(timeTime[0],BTime[0],s=50,marker='*',color='black')
    
    if(i == 0):
        ax1_2 = ax1.twinx()
    ax1_2.plot(timeTime,UcTime,color='xkcd:azure',linestyle=lStyle[i])
    ax1_2.scatter(timeTime[0],UcTime[0],s=50,marker='*',color='black')
    ax1_2.set_ylabel('Calving Rate [m/yr]',color='xkcd:azure')
    ax1_2.tick_params(axis='y',labelcolor='xkcd:azure')
    ax1.set_ylabel('Avg Melt Rate [m/day]',color='xkcd:red')
    ax1.tick_params(axis='y',labelcolor='xkcd:red')
    ax1.set_xlabel('Time [days]')
    ax1.grid(alpha=.5)
    ax1.set_title('Melt & Uc Forcing')
    if(args.labels != None):
        ax1.legend()

    ax2.plot(timeTime,pressTime,color='xkcd:sunflower',linestyle=lStyle[i])
    ax2.scatter(timeTime[0],pressTime[0],s=50,marker='*',color='black')
    ax2.set_ylabel('Mélange Force/Width [$Nm^{-1}$]',color='xkcd:black')
    ax2.tick_params(axis='y',labelcolor='xkcd:black')
    ax2.grid(alpha=.25)
    # ax2.set_yscale('log')
    if(i == 0):
        ax2_2 = ax2.twinx()
    ax2_2.plot(timeTime,VTime,color='xkcd:indigo',linestyle=lStyle[i],label='volume')
    ax2_2.set_ylabel('Mélange H0 [m]',color='xkcd:indigo')
    ax2_2.tick_params(axis='y',labelcolor='xkcd:indigo')
    ax2_2.set_ylabel('Volume [m^3]')
    # ax2_2.set_yscale('log')
    # ax2_2.grid(alpha=.25,color='xkcd:mulberry')
    ax2.set_xlabel('Time [days]')
    ax2.set_title('')

    if(i == 0): #Plot blanks for legend 
        ax3.plot([],[],color='xkcd:sand',linestyle='-',label='In')   
        ax3.plot([],[],color='xkcd:apple',linestyle='-',label='Out')  
        ax3.plot([],[],color='xkcd:red',linestyle='-',label='Melt')  
        ax3.plot([],[],color='xkcd:sky blue',linestyle='-',label='Net',alpha=.5)
        # ax3.plot([],[],color='xkcd:indigo',linestyle='-',label='Volume')
    ax3.plot(timeTime,inFluxTime,color='xkcd:sand',linestyle=lStyle[i])   
    ax3.plot(timeTime,endFluxTime,color='xkcd:apple',linestyle=lStyle[i])  
    ax3.plot(timeTime,bFluxTime,color='xkcd:red',linestyle=lStyle[i])
    ax3.plot(timeTime,inFluxTime-endFluxTime-bFluxTime,color='xkcd:sky blue',linestyle=lStyle[i],alpha=.5)
    ax3.plot([timeTime[0],timeTime[-1]],[0,0],color = 'gray',alpha=.5,linestyle='--')    
    ax3.set_ylabel('Volume Flux [m^3/s]')
    ax3.set_xlabel('Time [days]')
    ax3.grid(alpha=.5)
    ax3.set_title('')
    ax3.legend(loc='upper center', bbox_to_anchor=(0.45, 1.15),ncol=4)

    if(i == 0): #Plot blanks for legend 
        ax4.plot([],[],color='xkcd:apple',linestyle=lStyle[i],label='Out Flux Fraction')
        ax4.plot([],[],color='xkcd:red',linestyle=lStyle[i],label='Melt Flux Fraction')
        ax4.plot([],[],color='xkcd:sky blue',linestyle=lStyle[i],label='Balance Fraction',alpha=.5)
    ax4.plot(timeTime,endFluxTime/inFluxTime,color='xkcd:apple',linestyle=lStyle[i])
    ax4.plot(timeTime,bFluxTime/inFluxTime,color='xkcd:red',linestyle=lStyle[i])
    ax4.plot(timeTime,(bFluxTime + endFluxTime)/inFluxTime,color='xkcd:sky blue',linestyle=lStyle[i],alpha=.5)

    ax4.plot([timeTime[0],timeTime[-1]],[0,0],color = 'black',alpha=0)
    ax4.plot([timeTime[0],timeTime[-1]],[1,1],color = 'gray',alpha=.5,linestyle='--')
    ax4.set_ylabel('Fraction of Influx Volume')
    ax4.set_xlabel('Time [days]')
    ax4.grid(alpha=.5)
    ax4.set_title('')
    # ax4.legend()

dirStr = ''
if(os.path.isdir('figs')):
    dirStr = 'figs/'
multiStr = ''
if(i > 0):
    multiStr = f'Multi{i+1}'
plt.savefig(f'{dirStr}MelangeMass{multiStr}.png',format='png',dpi=150)
if(args.silent == 0):
    plt.show()
plt.close()
