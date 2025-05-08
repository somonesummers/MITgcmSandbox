
import numpy as np
from matplotlib import pyplot as plt
import sys
import fileinput
import shutil
from MITgcmutils import mds
import os
# sys.path.append('/hdd/glaciome/models/glaciome1D')
# sys.path.insert(0, '')
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import constants, glaciome
import glob
import pickle
import time


n_pts = 21 # number of grid points
L = 6e3 # ice melange length
Ut = 6e3 # glacier terminus velocity [m/a]; treated as a constant
Uc = 6e3 # glacier calving rate [m/a]; treated as a constant
Ht = 600 # terminus thickness
dt = 0.02# 1/(n_pts-1)/10 # time step [a]; needs to be quite small for this to work
# H0_manual = 94 #posit H0 start value, set to None to use default (75 m)

# specifying fjord geometry
X_fjord = np.linspace(0e3,400e3,101)
Wt = 4800
W_fjord = Wt + 0/10000*X_fjord
# B = -0.3*365
#already did 0.3, so skipping here
B_list = np.array([0.5,0.7,1.0,2.0]) * -1.0*365.0

for B in B_list:
	# first run to steady state
	data = glaciome(n_pts, dt, L, Ut, Uc, Ht, B, X_fjord, W_fjord, H0 = None)
	# files = sorted(glob.glob('couplingResults/MITgcmRun_[0-9][0-9][0-9][0-9][0-9].pickle'))
	# with open(files[-1], 'rb') as file:
	#     data = pickle.load(file)
	#     file.close()
	print(data)

	data.steadystate()
	print(data)
	data.save('SteadyUc%05iB%03iWt%05i.pickle' %(int(Uc),int(abs(B/365*100)),int(Wt)))
	os.system('mkdir Uc%05iB%03iWt%05i' %(int(Uc),int(abs(B/365*100)),int(Wt)))
	os.system('mv -v tempfile* Uc%05iB%03iWt%05i' %(int(Uc),int(abs(B/365*100)),int(Wt)))
