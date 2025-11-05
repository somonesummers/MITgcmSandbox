#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import sys
from scipy.integrate import simpson
import os
# sys.path.append('/hdd/glaciome/models/glaciome1D')
# sys.path.insert(0, '')
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import constants, glaciome

import glob
import pickle
import argparse

def makeColors(seedColor, n):
    color1 = 'xkcd:light ' + seedColor
    color2 = 'xkcd:dark ' + seedColor
    rlin = np.linspace(mcolors.to_rgb(color1)[0],mcolors.to_rgb(color2)[0],n)
    glin = np.linspace(mcolors.to_rgb(color1)[1],mcolors.to_rgb(color2)[1],n)
    blin = np.linspace(mcolors.to_rgb(color1)[2],mcolors.to_rgb(color2)[2],n)
    colors = np.stack((rlin,glin,blin), axis=0)
    return colors

# Take input options. Some defaults are set here, so be aware
parser = argparse.ArgumentParser(description='Plot options for melange profiles')
parser.add_argument('-f','--files', nargs=1, default=['couplingResults/MITgcmRun_'],
                    help='file string used to name pickle files FILES00000.pickle [default = couplingResults/MITgcm_]')
parser.add_argument('-u','--Uc', nargs=1, default=[0],type=int,
                    help='Plot Uc instead of side view [defaut = 0]')
parser.add_argument('-n','--NumView', nargs=1, default=[100], type=int,
                    help='How many files to look back for cascade plots [defaut = 100]')
parser.add_argument('-s','--silent', action='count', default=0,
                    help='Option to silence showing of plots')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
args = parser.parse_args()

# print(args)
fileStr = args.files[0]
files = sorted(glob.glob('%s*.pickle'%fileStr))

constant = constants()


# files = files[0:2]
# file = files[0]

# fig, axes = plt.subplots(2, 3, figsize=(12, 6), layout="constrained")
# ax1 = axes[0,0]
# ax2 = axes[0,1]
# ax3 = axes[1,0]
# ax4 = axes[1,1]
# ax5 = axes[0,2]
# ax6 = axes[1,2]

seedColor=['green','blue','violet','red']
# seedColor=['xkcd:dandelion','xkcd:rose','xkcd:lavender','xkcd:apple','xkcd:periwinkle','xkcd:seafoam','xkcd:umber']


jShift = 0
if(args.timeRange == None):
    toIterate = np.arange(len(files))
else:
    if(args.timeRange[1] == 0):
        toIterate = np.arange(len(files[args.timeRange[0]:]))
    else:
        toIterate = np.arange(len(files[args.timeRange[0]:args.timeRange[1]]))
    jShift = args.timeRange[0]

n = len(toIterate)
subSample = max([int(np.ceil(n / 10)),1])

H0Time = np.zeros(np.shape(toIterate))
VTime = np.zeros(np.shape(toIterate))
UcTime = np.zeros(np.shape(toIterate))
lengthTime = np.zeros(np.shape(toIterate))
timeTime = np.zeros(np.shape(toIterate))
iterationNumber = np.zeros(np.shape(toIterate))
meltTime = np.zeros([n,20])
for j in toIterate:
    file = files[j + jShift]
    with open(files[j + jShift], 'rb') as fileName:
        data = pickle.load(fileName)
        fileName.close()
    it = file.replace('./', '').replace('.pickle', '').replace(fileStr,'')
    H0Time[j]=data.H0
    X_ = np.concatenate(([data.X[0]], data.X_, [data.X[-1]]))
    H = np.concatenate(([data.H0], data.H, [data.HL]))
    W = np.concatenate(([data.W0], data.W, [data.WL]))
    VTime[j] = simpson(H*W, x=X_)*1e-9
    lengthTime[j]=data.X[-1] - data.X[0]
    UcTime[j] = data.Uc
    meltTime[j,:] = -1*data.B / 365.25
    timeTime[j] = data.t * 365.25
    iterationNumber[j]=int(it)

# seasonShift = 0
# colorList = ['xkcd:light blue','xkcd:apple','xkcd:tomato','xkcd:squash']
# plt.plot([],[],color=colorList[0],label='1')
# plt.plot([],[],color=colorList[1],label='2')
# plt.plot([],[],color=colorList[2],label='3')
# plt.plot([],[],color=colorList[3],label='4')
# for i in range(int(n/200)-1):
#     plt.plot(np.mean(meltTime[i*200+seasonShift:(i+1)*200+seasonShift,:],axis=0),color=colorList[i%4])
# plt.legend()
# plt.savefig(f'meltOverYears.png',format='png',dpi=200)
# plt.show()
# plt.close()

# seasonShift = 43
# colorList = ['xkcd:light blue','xkcd:apple','xkcd:tomato','xkcd:squash']
# plt.plot([],[],color=colorList[0],label='Winter')
# plt.plot([],[],color=colorList[1],label='Spring')
# plt.plot([],[],color=colorList[2],label='Summer')
# plt.plot([],[],color=colorList[3],label='Autumn')
# for i in range(int(n/91)-1):
#     plt.plot(np.mean(meltTime[i*91+seasonShift:(i+1)*91+seasonShift,:],axis=0),color=colorList[i%4])
# plt.legend()
# plt.savefig(f'meltOverSeasons.png',format='png',dpi=200)
# plt.show()
# plt.close()

colorList = ['xkcd:dark blue','xkcd:light blue',
                'xkcd:bright orange','xkcd:tomato',
                'xkcd:faded pink']
seasonShift = 21
fig, axes = plt.subplots(2, 3, figsize=(12, 6), layout="constrained")
axes = axes.reshape([6,1])[:,0]
axes[0].plot([],[],color=colorList[0],label='Early Winter')
axes[0].plot([],[],color=colorList[1],label='Late Winter')
axes[0].plot([],[],color=colorList[2],label='Plume Up')
axes[0].plot([],[],color=colorList[3],label='Plume')
axes[0].plot([],[],color=colorList[4],label='Plume down')
for i in range(int(n/73)-3):
    axes[0].plot(np.linspace(0,1,20),np.mean(meltTime[i*73+seasonShift:(i+1)*73+seasonShift,:],axis=0),color=colorList[i%5])
    axes[i%5 + 1].plot(np.linspace(0,1,20),np.transpose(meltTime[i*73+seasonShift:(i+1)*73+seasonShift,:]),color=colorList[i%5],alpha=.05)
for i in range(6):
    axes[i].set_ylim([0,1.6])
    if( i % 3 == 0):
        axes[i].set_ylabel('Net Melt Rate [m/day]')
    if(i > 2):
        axes[i].set_xlabel('Fraction along mélange [ ]')
axes[0].legend()
plt.savefig(f'figs/meltOverMonths.png',format='png',dpi=200)
plt.show()
plt.close()

# print("Making melange gif")
# os.system('magick -delay %f figs/advectBergs*.png -colors 256 -depth 256 figs/advectBergs.gif' %(500/n))

