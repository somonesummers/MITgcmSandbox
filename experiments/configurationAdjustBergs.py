#!/usr/bin/env python
# coding: utf-8

import fileinput
import shutil
import sys
import os
import glob

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

def changeDiagnosticFreq(oldFreq, newFreq):
    for i in range(4):
        replaceAll('input/data.diagnostics','frequency(%i)=%i,' %(i,oldFreq), 'frequency(%i)=%i,' %(i,newFreq))
    replaceAll('input/data.diagnostics','frequency(%i)=%i,' %(4,-oldFreq), 'frequency(%i)=%i,' %(4,-newFreq))

secsInADay = 24*3600

#Load old data file for editing
dt = 0.0   
nInterOld = 0
endTimeOld = 20*secsInADay
pickupIterationOld = 10*secsInADay
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
        elif "nIter0=" in line:
            nInterOld = float(line[8:-2])
        elif "endTime=" in line:
            endTimeOld = float(line[9:-2])
        elif "pChkptFreq=" in line:
            pickupIterationOld = float(line[12:-2])
print('dt is loaded as', dt)
print('nIterOld is loaded as', nInterOld)
print('endTimeOld is loaded as', endTimeOld)
print('pickupIterationOld is loaded as', pickupIterationOld)

# Standard spin-up times

endTimeNew = endTimeOld + 4*secsInADay
nInterNew = endTimeOld/dt
pickupIterationNew = .5*secsInADay

#set new iter0 value, adjust pickupsaving to 12 hours, extend runtime of sim
replaceAll('input/data','nIter0=%i,' % nInterOld,'nIter0=%i,' % nInterNew)
replaceAll('input/data','pChkptFreq=%.8e,' % pickupIterationOld, 'pChkptFreq=%.8e,' % pickupIterationNew)
replaceAll('input/data','endTime=%i' % endTimeOld, 'endTime=%i' % endTimeNew)

#Turn on Berg Diagnostics
replaceAll('input/data.diagnostics','# timePhase(2)', ' timePhase(2)')
replaceAll('input/data.diagnostics','# fields(1:3,3)', ' fields(1:3,3)')
replaceAll('input/data.diagnostics','# fileName(3)', ' fileName(3)')
replaceAll('input/data.diagnostics','# frequency(3)', ' frequency(3)')

#Make Diagnostics every 4 hours
changeDiagnosticFreq(.5*secsInADay,.25*secsInADay)

#This shouldn't be necessary anymore
#os.rename('input/data.iceberg_old','input/data.iceberg')

#Turn on ICEBERG package
replaceAll('input/data.pkg','useICEBERG=.FALSE.,', 'useICEBERG=.TRUE.,')

#bump results into holding folder
os.rename('results','results_init')
os.makedirs('results')

#Shuffle pickup file into input directory
dest_dir = 'input/'
for file in glob.glob('results_init/pickup.%010i.*' % nInterNew):
    print('Copying over pickup file:',file)
    shutil.copy(file, dest_dir)

print('Done prepping for Berg run ⁺₊❅.٩(^ᗜ^ )و')




