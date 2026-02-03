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
parser.add_argument('-f','--files', nargs=1, default=['couplingResults/MITgcm_'],
                    help='file string used to name pickle files FILES00000.pickle [default = couplingResults/MITgcm_]')
parser.add_argument('-u','--Uc', nargs=1, default=[0],type=int,
                    help='Plot Uc instead of side view [defaut = 0]')
parser.add_argument('-n','--numFrames', nargs=1, default=[60],type=int,
                    help='Number of Frames [defaut = 60]')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='Option to quickly make just last frame, qq shows it')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
args = parser.parse_args()
# print(args)
fileStr = args.files[0]
files = sorted(glob.glob('%s*.pickle'%fileStr))

constant = constants()


# files = files[0:2]
# file = files[0]



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
toIterate = toIterate + jShift
n = len(toIterate)
subSample = max([int(np.ceil(n / args.numFrames[0])),1])

toIterate = toIterate[::subSample]
n = len(toIterate)
shiftIndex = 0
print(n)
print(toIterate)

with open(files[0], 'rb') as file:
        print(f'initial file is {file}')
        data0 = pickle.load(file)
        file.close()


print('\tSubsampling at %i' %subSample)
# print(shiftIndex)
if(args.quick > 0):
    toIterate = [toIterate[-1]]
    subSample = 1
for j in toIterate[shiftIndex:]:
    fig, axes = plt.subplots(2, 2, figsize=(8, 6), layout="constrained")
    ax1 = axes[0,0]
    ax2 = axes[0,1]
    ax3 = axes[1,0]
    ax4 = axes[1,1]
    file = files[j]
    name = file.replace('./', '').replace('.pickle', '').replace(fileStr,'')
    # print(file)
    linestyle = '-'
    with open(files[j], 'rb') as file:
        data = pickle.load(file)
        file.close()
    X = data.X
    X_ = np.concatenate(([data.X[0]],data.X_,[data.X[-1]]))
    U = data.U
    H = np.concatenate(([data.H0],data.H,[1.5*data.H[-1]-0.5*data.H[-2]]))
    B = np.concatenate((data.B,[1.5*data.B[-1]-0.5*data.B[-2]])) * -1 #flip for plotting
    gg = np.concatenate(([1.5*data.gg[0]-0.5*data.gg[1]],data.gg,[1.5*data.gg[-1]-0.5*data.gg[-2]]))
    muW = data.muW# np.concatenate(([3*data.muW[0]-3*data.muW[1]+data.muW[2]],data.muW,[3*data.muW[-1]-3*data.muW[-2]+data.muW[-3]]))
    X = X-X[0]
    X_ = X_-X_[0]

    ax1.plot(X*1e-3,(U+data.Ut-data.Uc)/constant.daysYear,marker='o',color=seedColor[0],alpha=1,linestyle=linestyle,label=name)
    # ax1.plot(X*1e-3,(U+data.Ut-data.Uc)/constant.daysYear,marker='o',color=seedColor[j],linestyle=linestyle,label=names[j])
    ax1.set_xlabel('Distance Along Mélange [km]')
    ax1.set_ylabel('Speed [m/day]')
    ax1.grid(alpha=.5)
    # ax1.legend()

    if(args.Uc[0] == 0):
        ax2.plot([-2,20],[0,0],color='xkcd:ocean blue',linestyle='--',linewidth=0.5)
        ax2.plot(np.append(X_,X_[::-1])*1e-3,np.append(-constant.rho/constant.rho_w*H,(1-constant.rho/constant.rho_w)*H[::-1]),
            marker='o',color=seedColor[1],alpha=1,linestyle=linestyle,label=name)
        ax2.set_xlabel('Distance Along Mélange [km]')
        ax2.set_ylabel('Elevation [m]')
        # ax2.legend()
        ax2.set_xlim(ax1.get_xlim())
        ax2.grid(alpha=.5)

    ax3.plot(X_*1e-3,gg,marker='o',color=seedColor[2],alpha=1,linestyle=linestyle,label=name)
    ax3.set_xlabel('Distance Along Mélange [km]')
    ax3.set_ylabel('$g^{\\prime}$')
    ax3.grid(alpha=.5)

    ax4.plot(X_[1:-1]*1e-3,B[:-1]/constant.daysYear,marker='o',color=seedColor[3],alpha=1,linestyle=linestyle,label=name)
    ax4.plot([X_[1]*1e-3,X_[-1]*1e-3],[np.mean(B[:-1]/constant.daysYear),np.mean(B[:-1]/constant.daysYear)],color=seedColor[3],alpha=.25,linestyle='--',label=name)
    ax4.set_xlabel('Distance Along Mélange [km]')
    ax4.set_ylabel('Meltrate B [m/day]')     
    # ax4.legend()
    ax4.grid(alpha=.5)
    plt.suptitle(f'Mélange at day {(data.t- data0.t)*365.0:0.1f}')

    #Set limits of plots
    ax1.set_xlim([0,25])
    ax2.set_xlim([0,25])
    ax3.set_xlim([0,25])
    # ax4.set_xlim([0,25])

    ax1.set_ylim([50, 450])
    ax2.set_ylim(-200, 20)
    ax3.set_ylim([2,20])
    ax4.set_ylim([0,1.6])


    dirStr = ''
    timeStr = ''
    if(args.timeRange != None):
            timeStr = f'_T_{args.timeRange[0]}_{args.timeRange[1]}'
    if(os.path.isdir('figs')):
        dirStr = 'figs/'
    plt.savefig('%smelangeTemp%05i.png' %(dirStr,j),format='png',dpi=150)
    if(args.quick > 1):
        plt.show()
    plt.close()

if(args.quick == 0):
    os.system(f'magick -delay {500/n} {dirStr}melangeTemp*.png -colors 256 -depth 256 {dirStr}melange{timeStr}.gif')
    os.system(f'rm {dirStr}melangeTemp*.png')

