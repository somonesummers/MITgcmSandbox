#!/usr/bin/env python3

# Paul Summers November 2025
# Save SGD to file


import numpy as np
import sys
import argparse
import os
import fileinput
import pandas as pd
import xarray
from scipy import interpolate
from MITgcmutils import mds

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

x = mds.rdmds("results/XC") 
x = x - x[0,1]
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))

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

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)



dy = 2*y[0,0]
for line in fileinput.input(f'input/data'):
    if "ExternForcingPeriod=" in line:
        period = float(line[21:-2])
    elif "ExternForcingCycle=" in line:
        cycle = float(line[20:-2])
nt = int(cycle/period)
plumeMask = np.fromfile(f'input/plumeMask.bin', dtype='>f8') #(2 is line, 3 is semi-cone)
rad = np.fromfile(f'input/runoffRad.bin', dtype='>f8')
rad = rad.reshape((nt,len(plumeMask)))
vel = np.fromfile(f'input/runoffVel.bin', dtype='>f8')
vel = vel.reshape((nt,len(plumeMask)))

plumeType = np.max(plumeMask)
cleanRad = rad[:,plumeMask == plumeType]
cleanVel = vel[:,plumeMask == plumeType]
seasonTime = np.linspace(0,cycle/86400,nt)
if(plumeType == 2):
    runoff = dy*cleanRad[:,0]*cleanVel[:,0]
else: #semi-rad
    runoff = np.pi*cleanRad**2*cleanVel

f = interpolate.interp1d(seasonTime, runoff,fill_value='0')

time = np.arange(startStep,maxStep,sizeStep) * dt/86400

print(time)

time_coords = pd.to_timedelta(time, unit="D")
SGD = f(time)

# 2. Create the DataArray
da = xarray.DataArray(SGD, coords=[time_coords], dims=["time"])

# 3. Convert to a Dataset (optional, but standard practice for netCDFs)
ds = da.to_dataset(name="SGD over Time")
print(ds)
ds.to_netcdf('SGD.nc')







