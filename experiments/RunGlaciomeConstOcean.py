
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
import argparse

parser = argparse.ArgumentParser(description='Plot options for melange changed over time')
parser.add_argument('-f','--file', nargs='?', default=[None],
                    help='file string used to name initial pickle file FILES.pickle [default = None]')
parser.add_argument('-t','--timeRange', nargs='?', type=int, default = 1000,
                    help='Specification of days to run for [default = 1000]')
args = parser.parse_args()

def sysPrint(stringIn,tag): 
        print(stringIn)
        os.system('echo ""' + stringIn + '"" >> SOS_out' + tag + '.txt')


# # Load directly
files = sorted(glob.glob(args.file))
# files = sorted(glob.glob('SteadyUc06000B100Wt04800.pickle'))

with open(files[-1], 'rb') as file:
    data = pickle.load(file)
    file.close()
strTemp = f'loaded {file.name}'
sysPrint(strTemp,args.file.split('.')[0])

for i in range(args.timeRange):
	# first run to steady state
	data.prognostic()
	strTemp = f'iter {i}, H0: {data.H0} m, L: {data.L/1000} km'
	sysPrint(strTemp,args.file.split('.')[0])
	data.save(f'SOS_{args.file.split('.')[0]}_{i:05d}.pickle')

sysPrint('done!',args.file.split('.')[0])
