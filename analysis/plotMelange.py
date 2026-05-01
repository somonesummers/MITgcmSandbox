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
parser.add_argument('-ex','--extended', action='count', default=0,
                    help='Option to include the collapse tragectory')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
args = parser.parse_args()

# print(args)
fileStr = args.files[0]
files = sorted(glob.glob('%s*.pickle'%fileStr))
collapseIt = len(files) -1
files+= sorted(glob.glob('%s*.pickle'%(fileStr.replace('sweep','zTemp0'))))

constant = constants()


# files = files[0:2]
# file = files[0]

fig, axes = plt.subplots(2, 3, figsize=(12, 6), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax6 = axes[1,1]
ax5 = axes[0,2]
ax4 = axes[1,2]

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
    file = files[j + jShift]
    name = file.replace('./', '').replace('.pickle', '').replace(fileStr,'')
    print(file)
    linestyle = '-'
    with open(files[j + jShift], 'rb') as file:
        data = pickle.load(file)
        file.close()
    X = data.X
    X_ = np.concatenate(([data.X[0]],data.X_,[data.X[-1]]))
    U = data.U
    H = np.concatenate(([data.H0],data.H,[1.5*data.H[-1]-0.5*data.H[-2]]))
    B = np.concatenate((data.B,[1.5*data.B[-1]-0.5*data.B[-2]])) * -1 #flip for plotting
    gg = np.concatenate(([1.5*data.gg[0]-0.5*data.gg[1]],data.gg,[1.5*data.gg[-1]-0.5*data.gg[-2]]))
    g_loc = np.concatenate(([1.5*data.g_loc[0]-0.5*data.g_loc[1]],data.g_loc,[1.5*data.g_loc[-1]-0.5*data.g_loc[-2]]))
    # print(np.sqrt((np.diff(data.U)/np.diff(data.X))**2))
    # mu = np.sqrt((np.diff(U)/data.dx)**2)/(data.L*data.gg)
    mu = (np.sqrt((np.diff(data.U)/data.dx)**2)/data.L+data.param.deps)/(data.gg/data.param.Uscale*data.param.Lscale)
    muW = data.muW# np.concatenate(([3*data.muW[0]-3*data.muW[1]+data.muW[2]],data.muW,[3*data.muW[-1]-3*data.muW[-2]+data.muW[-3]]))
    X = X-X[0]
    X_ = X_-X_[0]

    colors = makeColors(seedColor[0],n)
    alphaList = np.arange(.1,.5,.4/n)
    alphaList[-subSample] = 1
    ax1.plot(X*1e-3,(U+data.Ut-data.Uc)/constant.daysYear,marker='o',color=colors[:,j-shiftIndex],alpha=alphaList[j-shiftIndex],
        linestyle=linestyle,label=name)
    # ax1.plot(X*1e-3,(U+data.Ut-data.Uc)/constant.daysYear,marker='o',color=seedColor[j],linestyle=linestyle,label=names[j])
    ax1.set_xlabel('Distance Along Mélange [km]')
    ax1.set_ylabel('Speed [m/day]')
    # ax1.set_title(f'Last {n} steps')
    ax1.grid(alpha=.5)
    # if j == shiftIndex:
    #     ax1.legend()

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
    ax3.set_ylabel('$g^{\\prime}$ [1/yr]')
    ax3.grid(alpha=.5)
    # ax3.legend()

    # ax4.plot(X*1e-3,muW,color='xkcd:lavender',linestyle=linestyle,label=names[j])
    # ax4.set_xlabel('Distance Along Fjord [km]')
    # ax4.set_ylabel('$\\mu_w$')  

    colors = makeColors(seedColor[3],n)   
    ax4.plot(X_[1:-1]*1e-3,B[:-1]/constant.daysYear,marker='o',color=colors[:,j-shiftIndex],alpha=alphaList[j-shiftIndex],linestyle=linestyle,label=name)
    ax4.set_xlabel('Distance Along Mélange [km]')
    ax4.set_ylabel('Meltrate B [m/day]')
    # ax4.legend()
    ax4.grid(alpha=.5)

    colors = makeColors('orange',n)
    colors2 = makeColors('plum',n)
    ax5.plot(X*1e-3,muW,marker='o',color=colors[:,j-shiftIndex],alpha=alphaList[j-shiftIndex],linestyle=linestyle,label=name)
    ax5.set_xlabel('Distance Along Mélange [km]')
    ax5.set_ylabel('$\\mu_w$ []',color='xkcd:orange')
    ax5.tick_params(axis='y',labelcolor='xkcd:orange')
    if( j == shiftIndex):
        ax5_2 = ax5.twinx()
    ax5_2.plot(data.X_*1e-3,mu,marker='o',color=colors2[:,j-shiftIndex],alpha=alphaList[j-shiftIndex],linestyle=linestyle,label=name)
    x_lims = ax5.get_xlim()
    ax5_2.plot(ax5.get_xlim(),[.3,.3],color='xkcd:plum',linestyle='--',alpha =.15)
    ax5_2.set_ylabel('$\\mu$ []',color='xkcd:plum')
    ax5_2.tick_params(axis='y',labelcolor='xkcd:plum')
    ax5.set_xlim(x_lims)
    ax5.grid(alpha=.5)

    colors = makeColors('periwinkle',n)
    ax6.plot(X_*1e-3,g_loc,marker='o',color=colors[:,j-shiftIndex],alpha=alphaList[j-shiftIndex],linestyle=linestyle,label=name)
    ax6.set_xlabel('Distance Along Mélange [km]')
    ax6.set_ylabel('g_loc [1/yr]')
    ax6.grid(alpha=.5)

figLabels = ['(a)','(b)','(c)','(d)','(e)','(f)']
for label,ax in zip(figLabels,axes.reshape(6,1)):
    ax=ax[0]
    ax.text(
        ax.get_xlim()[0], ax.get_ylim()[1], label,
        fontsize='x-large', va='bottom',ha='right', fontfamily='sans serif')



strTemp = 'View'
if(args.extended == 1):
    strTemp = 'ViewExtended'
dirStr = ''
if(os.path.isdir('figs')):
    dirStr = 'figs/'
plt.savefig('%smelange%s.png' %(dirStr,strTemp),format='png',dpi=150)
if(args.silent == 0):
    plt.show()
plt.close()

# print("Making melange gif")
# os.system('magick -delay %f figs/advectBergs*.png -colors 256 -depth 256 figs/advectBergs.gif' %(500/n))

