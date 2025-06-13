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
parser.add_argument('-n','--NumView', nargs=1, default=[100], type=int,
                    help='How many files to look back for cascade plots [defaut = 100]')
parser.add_argument('-s','--silent', action='count', default=0,
                    help='Option to silence showing of plots')
args = parser.parse_args()

# print(args)
fileStr = args.files[0]
files = sorted(glob.glob('%s*.pickle'%fileStr))

constant = constants()


# files = files[0:2]
# file = files[0]

fig, axes = plt.subplots(2, 3, figsize=(12, 6), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]
ax5 = axes[0,2]
ax6 = axes[1,2]

seedColor=['green','blue','violet','red']
# seedColor=['xkcd:dandelion','xkcd:rose','xkcd:lavender','xkcd:apple','xkcd:periwinkle','xkcd:seafoam','xkcd:umber']
n = len(files)
subSample = max([int(np.ceil(n / 10)),1])

toIterate = np.arange(n)
H0Time = np.zeros(np.shape(toIterate))
VTime = np.zeros(np.shape(toIterate))
UcTime = np.zeros(np.shape(toIterate))
lengthTime = np.zeros(np.shape(toIterate))
timeTime = np.zeros(np.shape(toIterate))
iterationNumber = np.zeros(np.shape(toIterate))
for j in toIterate:
    file = files[j]
    with open(files[j], 'rb') as fileName:
        data = pickle.load(fileName)
        fileName.close()
    it = file.replace('./', '').replace('.pickle', '').replace(fileStr,'')
    H0Time[j]=data.H0
    X_ = np.concatenate(([data.X[0]], data.X_, [data.X[-1]]))
    H = np.concatenate(([data.H0], data.H, [data.HL]))
    W = np.concatenate(([data.W0], data.W, [data.WL]))
    VTime[j] = simpson(H*W, x=X_)*1e-9
    lengthTime[j]=data.X[-1]
    UcTime[j] = data.Uc
    timeTime[j] = data.t * 365
    iterationNumber[j]=int(it)

if(timeTime[0] != 0):
    print(f"time shifted by {timeTime[0]:.02f} days")
    timeTime = timeTime - timeTime[0] #remove any shift so that time starts at 0

if(args.Uc[0] == 1):
    ax2.plot(iterationNumber,UcTime,'-')
    sca=ax2.scatter(iterationNumber,UcTime,s=None,c=timeTime,cmap='viridis')
    # cbar=plt.colorbar(sca)
    # cbar.set_label('Iteration [Days]')
    ax2.set_ylabel('Calving Speed [m/yr]')
    ax2.set_xlabel('Iteration [ ]')
    ax2.grid(alpha=.5)


ax5.plot(iterationNumber,lengthTime,'-')
sca=ax5.scatter(iterationNumber,lengthTime,s=None,c=timeTime,cmap='viridis')
cbar=plt.colorbar(sca)
cbar.set_label('Time [Days]')
ax5.set_ylabel('Mélange length [m]')
ax5.set_xlabel('Time [days]')
ax5.grid(alpha=.5)

ax6.plot(H0Time,lengthTime,'-')
sca=ax6.scatter(H0Time,lengthTime,s=None,c=timeTime,cmap='viridis')
cbar=plt.colorbar(sca)
cbar.set_label('Time [Days]')
ax6.set_ylabel('Mélange Length [m]')
ax6.set_xlabel('Mélange H0 [m]')
ax6.grid(alpha=.5)


shiftIndex = 0
limit = args.NumView[0]
if(n > limit):
    shiftIndex = toIterate[-limit]
    # toIterate = toIterate[-100::] - toIterate[-100]
    subSample = round(limit /10) #fraction we pick
    n = limit
    print('\t== Displaying Only Final %i Steps in Color Plots==' %limit)
    # print(toIterate)
elif(subSample > 1):
    print('\tSubsampling at %i' %subSample)
# print(shiftIndex)
for j in toIterate[shiftIndex::subSample]:
    file = files[j]
    name = file.replace('./', '').replace('.pickle', '').replace(fileStr,'')
    print(file)
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

    colors = makeColors(seedColor[0],n)
    alphaList = np.arange(.1,.5,.4/n)
    alphaList[-subSample] = 1
    ax1.plot(X*1e-3,(U+data.Ut-data.Uc)/constant.daysYear,marker='o',color=colors[:,j-shiftIndex],alpha=alphaList[j-shiftIndex],linestyle=linestyle,label=name)
    # ax1.plot(X*1e-3,(U+data.Ut-data.Uc)/constant.daysYear,marker='o',color=seedColor[j],linestyle=linestyle,label=names[j])
    ax1.set_xlabel('Distance Along Mélange [km]')
    ax1.set_ylabel('Speed [m/day]')
    ax1.grid(alpha=.5)
    # ax1.legend()

    if(args.Uc[0] == 0):
        colors = makeColors(seedColor[1],n)
        ax2.plot([-2,20],[0,0],color='xkcd:ocean blue',linestyle='--',linewidth=0.5)
        ax2.plot(np.append(X_,X_[::-1])*1e-3,np.append(-constant.rho/constant.rho_w*H,(1-constant.rho/constant.rho_w)*H[::-1]),
            marker='o',color=colors[:,j-shiftIndex],alpha=alphaList[j-shiftIndex],linestyle=linestyle,label=name)
        ax2.set_xlabel('Distance Along Mélange [km]')
        ax2.set_ylabel('Elevation [m]')
        # ax2.legend()
        ax2.set_xlim(ax1.get_xlim())
        ax2.grid(alpha=.5)

    colors = makeColors(seedColor[2],n)
    ax3.plot(X_*1e-3,gg,marker='o',color=colors[:,j-shiftIndex],alpha=alphaList[j-shiftIndex],linestyle=linestyle,label=name)
    ax3.set_xlabel('Distance Along Mélange [km]')
    ax3.set_ylabel('$g^{\\prime}$')
    ax3.grid(alpha=.5)
    # ax3.legend()

    # ax4.plot(X*1e-3,muW,color='xkcd:lavender',linestyle=linestyle,label=names[j])
    # ax4.set_xlabel('Distance Along Fjord [km]')
    # ax4.set_ylabel('$\\mu_w$')  

    colors = makeColors(seedColor[3],n)   
    ax4.plot(X*1e-3,B/constant.daysYear,marker='o',color=colors[:,j-shiftIndex],alpha=alphaList[j-shiftIndex],linestyle=linestyle,label=name)
    ax4.set_xlabel('Distance Along Mélange [km]')
    ax4.set_ylabel('Meltrate B [m/day]')     
    # ax4.legend()
    ax4.grid(alpha=.5)

strTemp = 'View'
if(args.Uc[0] == 1):
    strTemp = 'ViewUc'
dirStr = ''
if(os.path.isdir('figs')):
    dirStr = 'figs/'
plt.savefig('%smelange%s.png' %(dirStr,strTemp),format='png',dpi=150)
if(args.silent == 0):
    plt.show()
plt.close()

# print("Making melange gif")
# os.system('magick -delay %f figs/advectBergs*.png -colors 256 -depth 256 figs/advectBergs.gif' %(500/n))

