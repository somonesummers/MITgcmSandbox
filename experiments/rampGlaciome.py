
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


# n_pts = 21 # number of grid points
# L = 6e3 # ice melange length
# Ut = 6e3 # glacier terminus velocity [m/a]; treated as a constant
# Uc = 6e3 # glacier calving rate [m/a]; treated as a constant
# Ht = 600 # terminus thickness
# dt = 0.02# 1/(n_pts-1)/10 # time step [a]; needs to be quite small for this to work
# # H0_manual = 94 #posit H0 start value, set to None to use default (75 m)

# # specifying fjord geometry
# X_fjord = np.linspace(0e3,400e3,101)
# Wt = 4800
# W_fjord = Wt + 0/10000*X_fjord
# # B = -0.3*365
# #already did 0.3, so skipping here
# B_list = np.array([0.5,0.7,1.0,2.0]) * -1.0*365.0

# # Load directly
files = sorted(glob.glob('SteadyUc06000B030Wt04800.pickle'))
# files = sorted(glob.glob('SteadyUc06000B100Wt04800.pickle'))

with open(files[-1], 'rb') as file:
    data = pickle.load(file)
    file.close()

Ramp_tmp = np.arange(6000,2000,-2)
Ramp_list = np.append(Ramp_tmp, Ramp_tmp[::-1])

data.t = 0

print('Running %i steps' %(len(Ramp_list)))
for i in range(len(Ramp_list)):
	# first run to steady state
	# data.B = -1.0 * B_list[i] * np.ones(len(data.H)) * 365 # melt is negative
	data.Uc = Ramp_list[i]
	data.dt = 1.0/365.0 #1 day
	# dL = 10
	# counter = 0
	# while dL > 2.0: # To ensure we go slow enough, we wait for dl/dt < 2 m/day ('steady state' is < 0.03 m/day)
	# while counter < 1: # Can also just take explict number of steps
	# 	oldL = data.L
	# 	data.prognostic(method='lm')
	# 	print('Length: %8.2f, H0: %7.3f Uc: %6.1f, t:%5.1f' %(data.L,data.H0,data.Uc,data.t*365))
	# 	counter += 1
	# 	dL = np.abs(oldL - data.L)
	# if(counter > 1):
	# 	print('\tCOUNTER: %i' %counter)
	data.prognostic(method='lm')  #OR just run it once
	print('Length: %8.2f m, H0: %7.3f m, Uc: %6.1f m/yr, t: %6.1f days' %(data.L,data.H0,data.Uc,data.t*365))

	if(i % 10 == 0):
		data.save('rampFile%05i.pickle' %(i))

