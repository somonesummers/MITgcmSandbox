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
parser.add_argument('-x','--xVariable', nargs=1, type=int, default=[1],
                    help='what to plot against 0: xVar , 1: Time [default]')
parser.add_argument('-s','--silent', action='count', default=0,
                    help='Option to silence showing of plots')
parser.add_argument('-ms','--muS', action='count', default=0,
                    help='Option to show muS on forcing plots')
parser.add_argument('-xt','--xTerminus', action='count', default=0,
                    help='Option to show terminus location on plots')
args = parser.parse_args()

# print(args)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
ax1 = axes[0]
ax2 = axes[1]

lStyle = ['-','--',':','-.']

if(args.silent > 0): #print status for script
    print(f"rampPlot: {args}")

for i in range(len(args.files)):
    fileStr = args.files[i]
    files = sorted(glob.glob('%s*.pickle'%fileStr))
    toIterate = np.arange(len(files))
    VTime = np.zeros(np.shape(toIterate))
    iterationNumber = np.zeros(np.shape(toIterate))
    bFluxTime = np.zeros(np.shape(toIterate))
    timeTime = np.zeros(np.shape(toIterate))
    endFluxTime = np.zeros(np.shape(toIterate))
    inFluxTime = np.zeros(np.shape(toIterate))
    for j in toIterate:
        file = files[j]
        with open(files[j], 'rb') as fileName:
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
        timeTime[j] = data.t*365.0 

        iterationNumber[j]=int(it)
    timeTime = timeTime - timeTime[0] # shift to start at t=0

    plt.suptitle('Mélange Mass Over Time')
    if(i == 0): #Plot blanks for legend 
        ax1.plot([],[],color='xkcd:apple',linestyle=lStyle[i],label='Out Flux Fraction')
        ax1.plot([],[],color='xkcd:red',linestyle=lStyle[i],label='Melt Flux Fraction')
        ax1.plot([],[],color='xkcd:sky blue',linestyle=lStyle[i],label='Balance Fraction',alpha=.5)
        ax1.plot([],[],color='xkcd:indigo',label='Volume')
    ax1.plot(timeTime,endFluxTime/inFluxTime,color='xkcd:apple',linestyle=lStyle[i])
    ax1.plot(timeTime,bFluxTime/inFluxTime,color='xkcd:red',linestyle=lStyle[i])
    ax1.plot(timeTime,(bFluxTime + endFluxTime)/inFluxTime,color='xkcd:sky blue',linestyle=lStyle[i],alpha=.5)

    ax1.plot([timeTime[0],timeTime[-1]],[0,0],color = 'black',alpha=0)
    ax1.plot([timeTime[0],timeTime[-1]],[1,1],color = 'gray',alpha=.5,linestyle='--')
    if(i == 0):
        ax1_1 = ax1.twinx() #only make new axis first time through
    ax1_1.plot(timeTime,VTime,color='xkcd:indigo',linestyle=lStyle[i],label='volume')
    ax1_1.set_ylabel('Mélange H0 [m]',color='xkcd:indigo')
    ax1_1.tick_params(axis='y',labelcolor='xkcd:indigo')
    ax1_1.set_ylabel('Volume [m^3]')
    
    ax1.set_ylabel('Fraction of Influx Volume')
    ax1.set_xlabel('Time [days]')
    ax1.grid(alpha=.5)
    ax1.set_title('')
    ax1.legend()

    if(i == 0): #Plot blanks for legend 
        ax2.plot([],[],color='xkcd:sand',linestyle='-',label='In Flux')   
        ax2.plot([],[],color='xkcd:apple',linestyle='-',label='Out Flux')  
        ax2.plot([],[],color='xkcd:red',linestyle='-',label='Melt Flux')  
        ax2.plot([],[],color='xkcd:sky blue',linestyle='-',label='Diff',alpha=.5)
    ax2.plot(timeTime,inFluxTime,color='xkcd:sand',linestyle=lStyle[i])   
    ax2.plot(timeTime,endFluxTime,color='xkcd:apple',linestyle=lStyle[i])  
    ax2.plot(timeTime,bFluxTime,color='xkcd:red',linestyle=lStyle[i])
    ax2.plot(timeTime,inFluxTime-endFluxTime-bFluxTime,color='xkcd:sky blue',linestyle=lStyle[i],alpha=.5)
    ax2.plot([timeTime[0],timeTime[-1]],[0,0],color = 'gray',alpha=.5,linestyle='--')    
    ax2.set_ylabel('Volume Change [m^3/s]')
    ax2.set_xlabel('Time [days]')
    ax2.grid(alpha=.5)
    ax2.set_title('')
    ax2.legend()

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
