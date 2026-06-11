from MITgcmutils import mds
import os
import sys
import fileinput
import numpy as np
import xmitgcm ## this is new
import argparse

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

i_s = list(range(5760,2108160,5760))
ds = xmitgcm.open_mdsdataset('results/',ignore_unknown_vars=False, delta_t = 15,iters=i_s,geometry='cartesian')

print(ds)

#this enables compression for all variables (lossless), levels are 1-9 fastest to most compressed
encoding = {var: {"zlib": True, "complevel": 3} for var in ds.data_vars}
# print(ds.variables)

ds.to_netcdf('xarray_data.nc',encoding=encoding) #this is slow, be patient

print('Done saving!')

