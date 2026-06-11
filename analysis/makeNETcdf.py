from MITgcmutils import mds
import os
import sys
import fileinput
import numpy as np
import xmitgcm ## this is new

print('Making xArray datafile now')

dt = 0.0
oldStartIter = 0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
        if "nIter0" in line:
            oldStartIter = int(line[8:-2])

## Cleaning first time step to global files. This sometimes had a hiccup in runs in early Oct 2025
print(f'dt loaded as {dt}')

# ds = xmitgcm.open_mdsdataset('results/',ignore_unknown_vars=False, delta_t = dt)

ds = xmitgcm.open_mdsdataset('results/',ignore_unknown_vars=False, delta_t = 15,iters='all')

print(ds)

# print(ds.variables)

ds.to_netcdf('xarray_data.nc', encoding={"time": {"dtype": "int64"}}) #this is slow, be patient

print('Done saving!')

