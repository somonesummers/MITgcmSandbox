from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import cmocean
import gsw
import utm

import geopandas as gpd


file2read = 'SF2015CTD.nc'
file_id = Dataset(file2read)
depth = np.asarray(file_id.variables['z'])
temp = np.asarray(file_id.variables['temp'])
temp_qc = file_id.variables['temp']
salt = np.asarray(file_id.variables['sal'])
temp[salt < 0] = np.nan
salt[salt < 0] = np.nan

np.save('2015_Fjord_T',temp[:,19])
np.save('2015_Fjord_Z',np.asarray(depth[:]))

file2read = f'Sermilik2015_XCTD.nc'
file_id = Dataset(file2read)
depth = file_id.variables['depth']
temp = np.asarray(file_id.variables['temp'])
salt = np.asarray(file_id.variables['sal'])

np.save('2015_Melange_T',temp[:,3])
np.save('2015_Melange_Z',np.asarray(depth[:]))

file2read = f'Sermilik2010m_XCTD.nc'
file_id = Dataset(file2read)
depth = file_id.variables['depth']
temp = np.asarray(file_id.variables['temp'])
salt = np.asarray(file_id.variables['sal'])
print(np.shape(temp))
np.save('2010mar_Melange_T',temp[:,3])
np.save('2010mar_Fjord_T',temp[:,5])
np.save('2010mar_Z',np.asarray(depth[:]))