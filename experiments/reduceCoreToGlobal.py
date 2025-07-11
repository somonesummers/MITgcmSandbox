#!/usr/bin/env python3
"""
Created on Wed Mar 26 2025

@author: psummers8
"""

import numpy as np
from matplotlib import pyplot as plt
import sys
import fileinput
import shutil
from MITgcmutils import mds
import os

def sysPrint(stringIn):
        print(stringIn)
        os.system('echo ""' + stringIn + '"" >> couplingResults/out.txt')

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

maxStep = 0
sizeStep = 1e10
startStep = 1e10

prefixes = ['BRGFlx','dynDiag','ptraceDiag','plumeDiag']

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

for endIter in np.arange(startStep, maxStep + 1, sizeStep):
    for k in range(len(prefixes)):
        sysPrint('\t\t == %s ==' %prefixes[k])
        dataTemp = mds.rdmds("results/%s"%(prefixes[k]), endIter)
        mds.wrmds('results/%s' %prefixes[k],dataTemp,itr=endIter, dataprec='float32')
        os.system('rm results/%s.%010i.0*.0*' %(prefixes[k],endIter))

prefixes = ['Depth','DXC','DXF','DXG','DXV','DYC','DYF','DYG','DYU','hFacC','hFacS','hFacW',
                    'maskInC','maskInS','maskInW','RAC','RAS','RAW','RAZ','XC','XG','YC','YG']
sysPrint('\tcondensing grid tile files to global files')
for k in range(len(prefixes)):
    sysPrint('\t\t== %s ==' %prefixes[k])
    dataTemp = mds.rdmds("results/%s"%(prefixes[k]))
    mds.wrmds('results/%s' %prefixes[k],dataTemp,dataprec='float32')
    os.system('rm results/%s.0*.0*' %(prefixes[k])) #picks out tile level files

sysPrint('Done!')