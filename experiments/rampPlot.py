#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import sys
from scipy.integrate import simpson
import argparse

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
parser.add_argument('-f','--files', nargs=1, default=['rampFile'],
                    help='file string used to name pickle files FILES00000.pickle [default = rampFile]')
parser.add_argument('-x','--xVariable', nargs=1, type=int, default=[0],
                    help='what to plot against 0: xVar [default], 1: Time')
parser.add_argument('-r','--rampVariable', nargs=1, type=int, default=[0],
                    help='what is ramping 0: MeltRate [default], 1: Uc')
args = parser.parse_args()

# print(args)
fileStr = args.files[0]
files = sorted(glob.glob('%s*.pickle'%fileStr))

# files = files[0:2]
# file = files[0]

fig, axes = plt.subplots(2, 2, figsize=(10, 6), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[1,0]
ax3 = axes[0,1]
ax4 = axes[1,1]

toIterate = np.arange(len(files))
H0Time = np.zeros(np.shape(toIterate))
VTime = np.zeros(np.shape(toIterate))
UcTime = np.zeros(np.shape(toIterate))
lengthTime = np.zeros(np.shape(toIterate))
iterationNumber = np.zeros(np.shape(toIterate))
pressTime = np.zeros(np.shape(toIterate))
BTime = np.zeros(np.shape(toIterate))
gTime = np.zeros(np.shape(toIterate))
spdTime = np.zeros(np.shape(toIterate))
timeTime = np.zeros(np.shape(toIterate))
for j in toIterate:
    file = files[j]
    with open(files[j], 'rb') as fileName:
        data = pickle.load(fileName)
        fileName.close()
    it = it = file.replace('./', '').replace('.pickle', '').replace(fileStr,'')
    H0Time[j]=data.H0
    X_ = np.concatenate(([data.X[0]], data.X_, [data.X[-1]]))
    H = np.concatenate(([data.H0], data.H, [data.HL]))
    W = np.concatenate(([data.W0], data.W, [data.WL]))
    VTime[j] = simpson(H*W, x=X_)*1e-9
    lengthTime[j]=data.X[-1]
    UcTime[j] = data.Uc
    BTime[j] = -1 * np.mean(data.B) /365.0 # we flip this for plotting purposes
    timeTime[j] = data.t*365.0 
    pressTime[j] = data.H0*data.pressure(data.H0)
    spdTime[j]=np.mean(data.U)
    gTime[j]=np.mean(data.gg)
    iterationNumber[j]=int(it)

# timeTime = (timeTime - timeTime[0])/5 #if pre-code fix for total time
if(args.rampVariable[0] == 0):
    forceTime = BTime
    forceLabel = 'Melt Rate [m/day]'
    forceString = 'Melt'
elif(args.rampVariable[0] == 1):
    forceTime = UcTime
    forceLabel = 'Calving Rate [m/yr]'
    forceString = 'Calving'
else:
    raise Exception('Invalid forcing option %i' %args.deltaVariable[0])


plt.suptitle('Mélange with Time Varying %s' %forceString)

ax1.plot(timeTime,forceTime,'-',color='xkcd:tomato')
ax1.scatter(timeTime[0],forceTime[0],s=50,marker='*',color='black')
# sca=ax1.scatter(BTime,H0Time,s=None,c=timeTime,cmap='cividis')
# cbar=plt.colorbar(sca)
# cbar.set_label('Iteration')
ax1.set_ylabel(forceLabel)
ax1.set_xlabel('Time [days]')
ax1.grid(alpha=.5)
ax1.set_title('%s Forcing' %forceString)

xString =''
if(args.xVariable[0] == 0):
    drivingVariable = forceTime
    drivingLabel = forceLabel
    xString = forceString
elif(args.xVariable[0] == 1):
    drivingVariable = timeTime
    drivingLabel = 'Time [days]'
    xString = 'Time'
else:
    raise Exception('Invalid X-axis option: %s' %args.xVariable[0])

ax2.plot(drivingVariable,lengthTime,'-',color='xkcd:blueberry')
ax2.scatter(drivingVariable[0],lengthTime[0],s=50,marker='*',color='black')
ax2.set_ylabel('Mélange length [m]',color='xkcd:blueberry')
ax2.tick_params(axis='y',labelcolor='xkcd:blueberry')
ax2.grid(alpha=.25,color='xkcd:blueberry')
# sca=ax2.scatter(BTime,lengthTime,s=None,c=timeTime,cmap='cividis')
# cbar=plt.colorbar(sca)
# cbar.set_label('Iteration')
ax2_2 = ax2.twinx()
ax2_2.plot(drivingVariable,H0Time,'-',color='xkcd:mulberry')
ax2_2.scatter(drivingVariable[0],H0Time[0],s=50,marker='*',color='black')
ax2_2.set_ylabel('Mélange H0 [m]',color='xkcd:mulberry')
ax2_2.tick_params(axis='y',labelcolor='xkcd:mulberry')
ax2_2.grid(alpha=.25,color='xkcd:mulberry')
ax2.set_xlabel(drivingLabel)
ax2.set_title('Size')

ax3.plot(drivingVariable,pressTime,'-',color='xkcd:sunflower')
ax3.scatter(drivingVariable[0],pressTime[0],s=50,marker='*',color='black')
# sca=ax3.scatter(BTime,pressTime,s=None,c=timeTime,cmap='cividis')
# cbar=plt.colorbar(sca)
# cbar.set_label('Time [days]')
ax3.set_ylabel('Mélange F/W [$Nm^{-1}$]')
ax3.set_yscale('log')
ax3.set_xlabel(drivingLabel)
ax3.grid(alpha=.5)
ax3.set_title('Buttressing Strength')

ax4.plot(drivingVariable,spdTime/365.0,'-',color='xkcd:apple')
ax4.scatter(drivingVariable[0],spdTime[0]/365.0,s=50,marker='*',color='black')
ax4.set_ylabel('Speed [m/day]',color='xkcd:apple')
ax4.tick_params(axis='y',labelcolor='xkcd:apple')
ax4.grid(alpha=.25,color='xkcd:apple')
ax4_2=ax4.twinx()
# sca=ax4.scatter(spdTime/365.0,gTime,s=20,c=timeTime,cmap='cividis')
# cbar=plt.colorbar(sca)
# cbar.set_label('Time [days]')
ax4_2.plot(drivingVariable,gTime,'-',color='xkcd:lilac')
ax4_2.scatter(drivingVariable[0],gTime[0],s=50,marker='*',color='black')
ax4_2.set_ylabel('$g^{\\prime}$    ',color='xkcd:lilac')
ax4_2.tick_params(axis='y',labelcolor='xkcd:lilac')
ax4_2.grid(alpha=.25,color='xkcd:lilac')
ax4.set_title('Speed and $g^{\\prime}$')
ax4.set_xlabel(drivingLabel)

plt.savefig('autoRamp%sResults%s.png' %(forceString,xString),format='png',dpi=150)
plt.show()
plt.close()
