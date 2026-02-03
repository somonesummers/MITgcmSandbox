#!/usr/bin/env python3

# Paul Summers Dec 2025
# Made for making GIF for AGU talk to match prettyPlot figure style

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
from matplotlib import colormaps
import sys
from scipy.integrate import simpson
import argparse
import os
import fileinput
# sys.path.append('/hdd/glaciome/models/glaciome1D')
# sys.path.insert(0, '')
sys.path.append('/Users/psummers8/Documents/glaciome1D')
from glaciome1D import constants, glaciome
from MITgcmutils import mds
import cmocean
import glob
import pickle
import gsw
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings

cleanPNGs = True
constant = constants()

# Take input options. Some defaults are set here, so be aware
parser = argparse.ArgumentParser(description='Plot options for melange changed over time')
parser.add_argument('-f','--files', nargs='*', default=['couplingResults/MITgcmRun_'],
                    help='file string used to name pickle files FILES00000.pickle [default = couplingResults/MITgcmRun_]')
parser.add_argument('-l','--labels', nargs='*', default=None,
                    help='file string label plots [default = None]')
parser.add_argument('-sgd','--sgd', nargs='*', type=float, default = None,
                    help='SGD, any value loads SGD from input, needs results RY')
parser.add_argument('--silent', action='count', default=1,
                    help='Option to silence showing of plots')
parser.add_argument('-n','--numFrames', nargs='?', type=int, default = 60,
                    help='optional specification of numFrames [default = 60]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = [94, 459],
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('-w','--windowMean', nargs='?', type=int, default=10,
                    help='window width for time averaging melt rates [default = 10]')
parser.add_argument('-ff','--fixForce', action='count', default=0,
                    help='Option to fix forcing left scale to +/- 5%%')
parser.add_argument('-ll','--legLoc', nargs='?', default=0,type=int,
                    help='legend location [default = 0]')
parser.add_argument('-dpi','--dpi', nargs='?', default=200,type=int,
                    help='dpi to print [default = 200]')
parser.add_argument('-y','--yCrossSection', nargs=1, type=float,default = 5000,
                    help='Y slice location [default 5000 m]')
parser.add_argument('-s','--shadow', action='count', default=1,
                    help='option of shadow for mélange [default (off)]')
parser.add_argument('-xl','--xlimit', nargs='?', type=int, default = 30,
                    help='optional max x [default = 30 km]')
parser.add_argument('-Y','--Years', nargs='*', type=int, default = [1,2],
                    help='years to plot, default = [1,2]')
args = parser.parse_args()

showDensity = False
showZero = False
showQuiver = True

# print(args)


lStyle = ['-',':','--','-.',(0, (3, 2, 1, 2, 1, 2)),(0, (2, 3, 1, 2, 1, 2, 1, 2,))]
meltColors = ['xkcd:tomato red','xkcd:rose','xkcd:wine','xkcd:grape purple','xkcd:grape','xkcd:violet']

if(args.silent > 0): #print status for script
    print(f"rampPlot: {args}")

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

maxStep = 0
sizeStep = 1e10
startStep = 1e10
for file in os.listdir('results'):
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
    maxStep = args.timeRange[1] * 86400 / dt

if((maxStep-startStep)/sizeStep > args.numFrames):   #if more than numFrames, downscale to be less than numFrames
    dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/args.numFrames)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

x = mds.rdmds("results/XC")/1e3 
x = x - x[0,1]
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))

## Do temp, U plotting first, only for first year
ySlice = np.argmin(np.abs(y[:,0] - args.yCrossSection))
print('cross section is y =', y[ySlice,0], 'index', ySlice)

for i in np.arange(startStep, maxStep + 1, sizeStep):
    plt.figure(figsize=(12, 4),layout="tight")
    ax5 = plt.gca()
    data = mds.rdmds("results/dynDiag", i)
    data[:,data[1,:,:,:] == 0] = np.nan #if salt = 0 its a wall, nan it out for averging later

    if(args.shadow > 0): #enable berg shadows here
        dataBergs = mds.rdmds("results/BRGFlx",i)
        if(dataBergs.shape[0] < 6):
            dataBergPlot = dataBergs[0,:,:,:].copy()
            dataBergPlot[dataBergPlot != 0] = .2
            dataBergPlot[data[0,:,:,:] < -1] = .2 #attemp to recapture freezing limit ice
            dataBergPlot[:,np.nansum(dataBergs[0,:,:,:],axis=0) == 0] = 0 #attemp to recapture freezing limit ice
        else:
            dataBergPlot = dataBergs[5,:,:,:]
        dataBergPlot[dataBergPlot == 0] = 1 #np.nan
    ## Temp
    lvl = np.linspace(-2, 4.0, 31)
    cm = 'RdYlBu_r'
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        plotData = np.squeeze(np.nanmean(data[0, :, :, :],axis=1))
    cp = ax5.contourf(
                    np.squeeze(x[ySlice,:]),
                    np.squeeze(z),
                    plotData,
                    lvl,
                    extend="both",
                    cmap=cm,
                    )
    if(True and showQuiver):
        zSkip = 1
        ySkip = 5
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            u = np.squeeze(np.nanmean(data[2, ::zSkip, :, ::ySkip],axis=1))
            w = np.squeeze(np.nanmean(data[3, ::zSkip, :, ::ySkip],axis=1))*10
        qv = ax5.quiver(
            x[ySlice,::ySkip],
            np.squeeze(z[::zSkip]),
            u,#/np.sqrt(u**2 + w**2 + 1e-12),
            w,#/np.sqrt(u**2 + w**2 + 1e-12),
            alpha=.4,
            #width = .004, #width of line
            scale = 1.2e-0
            )
        plt.quiverkey(qv,1.08,.95,.1,'10 cm/s')
    if(args.shadow > 0):
                # plt.plot(x[ySlice,:],-np.max(maxDepth,axis=0),color='gray',linestyle='dotted')
                cp2 = ax5.contourf(
                    x[ySlice,:],
                    np.squeeze(z),
                    np.squeeze(1-dataBergPlot[:, ySlice, :])*100,
                    [5,20,40,60,80],
                    extend="max",
                    alpha=.3,
                    cmap='cmo.gray_r')
    if(showDensity):
        salt = np.squeeze(data[1,:,ySlice,:])
        pressure = -1 * np.ones(salt.shape) * 1020 * 9.81 * np.repeat(np.expand_dims(z,1), salt.shape[1], axis=1) /10e3
        CT = gsw.CT_from_t(salt, data[0,:,ySlice,:], pressure)
        density = gsw.rho(salt, CT, 0) - 1000 #in-stu density less 1000
        densityLevels = np.linspace(22,28,31)
        cc = ax5.contour(
            np.squeeze(x[ySlice,:]),
            np.squeeze(z),
            np.squeeze(density),
            densityLevels,
            colors='black',
            linewidths=0.5,
            alpha=0.5
        )
        plt.clabel(cc, inline=3, fontsize=8)
    if(args.xlimit != None):
        ax5.set_xlim([0, args.xlimit])
    divider = make_axes_locatable(ax5)
    ax5_cb = divider.append_axes("right", size="2%", pad="15%")
    ax5_cb2 = divider.append_axes("right", size="2%", pad="10%")
    fig = ax5.get_figure()
    fig.add_axes(ax5_cb)
    fig.add_axes(ax5_cb2)
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")
    ax5.set_xlabel('Along fjord [km]')
    ax5.set_ylabel('Depth [m]')
    ax5.set_title(f"Fjord on day {(i-startStep)*dt/86400:4.0f}")
    cbar = plt.colorbar(cp,cax=ax5_cb,ticks=[-1, 0, 1,2,3,4], orientation="vertical",fraction=0.06,format='%.0f')
    cbar.set_label('Temperature [C]')
    cbar2 = plt.colorbar(cp2,cax=ax5_cb2, orientation="vertical",fraction=0.06,format='%.0f')
    cbar2.set_label('Ice fraction $\\varphi$ [%]')

    dirStr = ''
    if(os.path.isdir('figs')):
        dirStr = 'figs/'
    multiStr = ''
    timeString = ''
    if(args.timeRange != None):
        multiStr = f'{multiStr}_T_{args.timeRange[0]}_{args.timeRange[1]}'
    j = i/sizeStep + startStep
    plt.savefig('%sprettyQuiver%05d.png' %(dirStr,j),format='png',dpi=args.dpi)
    if(args.silent > 1):
        plt.show()
    plt.close()

os.system(f'magick -delay {800/((maxStep-startStep)/sizeStep)} figs/prettyQuiver*.png -colors 256 -depth 256 figs/prettyQuiver{multiStr}.gif')

#Clean up intermediate pngs
if(cleanPNGs):
    os.system('rm -f figs/prettyQuiver*.png')