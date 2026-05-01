from MITgcmutils import mds
import os
import sys
import fileinput
import numpy as np
import xmitgcm

# Clean tile level files
dt = 0.0
oldStartIter = 0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
        if "nIter0" in line:
            oldStartIter = int(line[8:-2])

## Cleaning first time step to global files. This sometimes had a hiccup in runs in early Oct 2025
print(f'dt loaded as {dt}')

prefixes = ['dynDiag','ptraceDiag','plumeDiag','presDiag','momDiag','heatDiag']
# prefixes = ['BRGFlx','dynDiag','ptraceDiag','plumeDiag','presDiag','momDiag','heatDiag']
# prefixes = ['presDiag','momDiag','heatDiag']

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

print(f'start {startStep}, step {sizeStep}, max {maxStep}')

for i in np.arange(startStep, maxStep + 1, sizeStep):
	print(f'condensing diagnostic tile files to global files for iteration {i}')
	# sysPrint('\tcondensing diagnostic tile files to global files iter:%i' %endIter)
	for k in range(len(prefixes)):
	    print('\t == %s ==' %prefixes[k])
	    dataTemp, it, meta = mds.rdmds("results/%s"%(prefixes[k]), i, returnmeta=True)
	    # print(meta)
	    if('fldlist' not in meta):
	    	print('Backfilling diagnostic metadata')
	    	if(prefixes[k] == 'BRGFlx'):
	    		meta['fldlist'] = ['BRGfwFlx','BRGhtFlx','BRGmltRt','BRG_TauX','BRG_TauY','BRGhFacC']
	    	elif(prefixes[k] == 'dynDiag'):
	    		meta['fldlist'] = ['THETA','SALT','UVEL','WVEL','VVEL']
	    	elif(prefixes[k] == 'ptraceDiag'):
	    		meta['fldlist'] = ['TRAC01','TRAC02']
	    	elif(prefixes[k] == 'plumeDiag'):
	    		meta['fldlist'] = ['icefrntW','icefrntT','icefrntS','icefrntA','icefrntR']
	    	elif(prefixes[k] == 'momDiag'):
	    		meta['fldlist'] = ['Um_Diss','Um_Advec','Um_Cori','Um_dPhiX','Wm_Diss','Wm_Advec']
	    	elif(prefixes[k] == 'heatDiag'):
	    		meta['fldlist'] = ['ADVx_TH','ADVy_TH','ADVr_TH','UTHMASS','VTHMASS','WTHMASS']
	    	elif(prefixes[k] == 'presDiag'):
	    		meta['fldlist'] = ['PHIHYD','PHI_NH']
	    mds.wrmds('results/%s' %prefixes[k],dataTemp,itr=i, dataprec='float32',fields=meta['fldlist'],)
	    os.system('rm results/%s.%010i.0*.0*' %(prefixes[k],i)) #picks out tile level files

# Clean up final files
prefixes = ['Depth','DXC','DXF','DXG','DXV','DYC','DYF','DYG','DYU','hFacC','hFacS','hFacW',
                    'maskInC','maskInS','maskInW','RAC','RAS','RAW','RAZ','XC','XG','YC','YG']
print('condensing grid tile files to global files')
for k in range(len(prefixes)):
    print('\t== %s ==' %prefixes[k])
    dataTemp = mds.rdmds("results/%s"%(prefixes[k]))
    mds.wrmds('results/%s' %prefixes[k],dataTemp,dataprec='float32')
    os.system('rm results/%s.0*.0*' %(prefixes[k])) #picks out tile level files

if(False):
	print('Making xArray datafile now')

	ds = xmitgcm.open_mdsdataset('results/',ignore_unknown_vars=False, delta_t = dt)

	print(ds)

	ds.to_netcdf('xarray_data.nc') #this is slow, be patient

print('Done saving!')

