#!/usr/bin/env python
# coding: utf-8

# In[2]:


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

#Load dt for this model run
dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            print(line)
            dt = float(line[8:-2])
print('dt is loaded as', dt)

# Standard spin-up times
endTimeOld = 24*secsInADay
endTimeNew = endTimeOld + 10*secsInADay

nInterOld = 20*secsInADay/dt
nInterNew = endTimeOld/dt

pickupIterationOld = .5*secsInADay
pickupIterationNew = 10*secsInADay

#set new iter0 value, adjust pickupsaving to 12 hours, extend runtime of sim
replaceAll('input/data','nIter0=%i,' % nInterOld,'nIter0=%i,' % nInterNew)
replaceAll('input/data','pChkptFreq=%.8e,' % pickupIterationOld, 'pChkptFreq=%.8e,' % pickupIterationNew)
replaceAll('input/data','endTime=%i' % endTimeOld, 'endTime=%i' % endTimeNew)

#Turn off Berg Diagnostics
replaceAll('input/data.diagnostics',' timePhase(2)', '# timePhase(2)')
replaceAll('input/data.diagnostics',' fields(1:3,3)', '# fields(1:3,3)')
replaceAll('input/data.diagnostics',' fileName(3)', '# fileName(3)')
replaceAll('input/data.diagnostics',' frequency(3)', '# frequency(3)')

#Make Diagnostics every hours
changeDiagnosticFreq(.25*secsInADay,.5*secsInADay)

#This shouldn't be necessary anymore
#os.rename('input/data.iceberg_old','input/data.iceberg')

#Turn off ICEBERG package
replaceAll('input/data.pkg','useICEBERG=.TRUE.,', 'useICEBERG=.FALSE.,')

#bump results into holding folder
os.rename('results','results_berg')
os.makedirs('results')

#Shuffle pickup file into input directory
dest_dir = 'input/'
for file in glob.glob('results_berg/pickup.%010i.*' % nInterNew):
    print('Copying over pickup file:'file)
    shutil.copy(file, dest_dir)

print('Done prepping for Post Berg run 𖤓°⋆.ೃ࿔･٩(^ᗜ^ )و')



