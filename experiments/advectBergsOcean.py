import os
import sys
import numpy as np
import fileinput
from MITgcmutils import mds
import matplotlib.pyplot as plt
import warnings
from scipy.interpolate import interp1d
import cmocean
import time
import random

start_time = time.time()

def write_bin(fname, data):
    # print(fname + " " + str(np.shape(data)))
    data.astype(">f8").tofile('input/'+fname)

## Controllers for which steps of advection process we do

advectBergs = True
clipBergs = True
trialOnly = False   
showPlots = False
addBergs = True
couplingTimeStep = 6*3600 #[s]

iceDensity = 917 # [kg/m^3]
oceanDensity = 1030 # [kg/m^3]

## Some rules on advecting/melting/sizing
minBergDepth = 5 # bergs not allowed to shrink/melt smaller than this [m]
maxBergDepth = 450 # bergs not allowed to shrink/melt smaller than this [m]
maxLambda = .90 # Cells are full at this lambda [hard limit]
## Bergs added to top off to icebergCoverLamba left of refreshGate
icebergRefreshingGate = 2
maxBergs = 500

## Load model wide parameters

sys.path.append('.')
from melangeModel import *

icebergPhi = np.min([.8,icebergCoverLambda]) # Scaling from GLACIOME1D
icebergCoverLambda = np.min([.9,icebergCoverLambda]) # How full filled cells are

print('\t\tLambda/Phi %.3f/%.3f' %(icebergCoverLambda,icebergPhi))
## Get most recent data from MITgcm

maxStep = 0
for file in os.listdir('results'):
    # print(file)
    if "BRGFlx.0" in file:
        words = file.split(".")
        # print(words[1])  
        if int(words[1]) > maxStep:
            maxStep = int(words[1])

## Iceberg melt data
dataMITgcm = mds.rdmds("results/BRGFlx", maxStep)
dataMITgcmOcean = mds.rdmds("results/dynDiag", maxStep)

## Geometry of MITgcm domain
x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = -1*np.squeeze(mds.rdmds("results/RC"))

ny = np.squeeze(np.shape(x)[0]) #x,y are loaded as [ny,nx] arrays from MITgcm
nx = np.squeeze(np.shape(x)[1])
nz = np.squeeze(np.shape(z)[0])
dz = np.load('input/dz.npy')
deltaX = x[0,0]*2
deltaY = y[0,0]*2

sum_z = np.cumsum(dz)

## This needs to know about z. Bergs cant be deeper than this
hardMaxDepth = np.max(np.abs(z)) - 5
## icebergs are destroyed when they move past here
icebergRightHandGate = int(nx*0.95)

## Now we load the iceberg geometry and mask files

## 2D Fields 
bathymetry  = np.fromfile('input/bathymetry.bin', dtype='>f8')
bathymetry  = bathymetry.reshape((ny,nx))
plumeMask  = np.fromfile('input/plumeMask.bin', dtype='>f8')
plumeMask  = plumeMask.reshape((ny,nx))
bergWidths = np.fromfile('input/icebergs_widths.bin', dtype='>f8')
bergWidths = bergWidths.reshape((maxBergs,ny,nx))
bergLength = np.fromfile('input/icebergs_length.bin', dtype='>f8')
bergLength = bergLength.reshape((maxBergs,ny,nx))
bergDepths = np.fromfile('input/icebergs_depths.bin', dtype='>f8')
bergDepths = bergDepths.reshape((maxBergs,ny,nx))
bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
bergMask = bergMask.reshape((ny,nx))
bergMaskNums = np.fromfile('input/bergMaskNums.bin', dtype='>f8')
bergMaskNums = bergMaskNums.reshape((ny,nx))
bergsPerCell = np.fromfile('input/numBergsPerCell.bin', dtype='>f8')
bergsPerCell = bergsPerCell.reshape((ny,nx)).astype(np.int64) 
driftMask = np.fromfile('input/driftMask.bin', dtype='>f8')
driftMask = driftMask.reshape((ny,nx))
barrierMask = np.fromfile('input/barrierMask.bin', dtype='>f8')
barrierMask = barrierMask.reshape((ny,nx))
meltMask = np.fromfile('input/meltMask.bin', dtype='>f8')
meltMask = meltMask.reshape((ny,nx))
## 3D Fields
totalBergArea = np.fromfile('input/totalBergArea.bin', dtype='>f8')
totalBergArea = totalBergArea.reshape((nz,ny,nx))
openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
openFrac = openFrac.reshape((nz,ny,nx))

# lam = np.nanmean((1-openFrac[0,1:-1,:]),axis=0)
# lam2 = np.mean(np.nansum(bergDepths[:,1:-1,:]*bergWidths[:,1:-1,:]*bergLength[:,1:-1,:],axis=0)/(deltaX*deltaY),axis=0)
# print('lambda',lam,lam2)

## Create fresh of geometry arrays, this lets us move bergs without double counting

bergWidthsNew   = np.zeros_like(bergWidths)
bergLengthNew   = np.zeros_like(bergWidths)
bergDepthsNew   = np.zeros_like(bergWidths)
bergsPerCellNew = np.zeros_like(bergsPerCell).astype(np.int64)

print('\t\tTotal bergs at start move: %i' %np.sum(bergsPerCell))

## Advect bergs. We go backwards (front to back) to ensure overfull icebergs can get priority to stay where they are
## If you want it to not really be random
# np.random.seed(3)
maxMoveX = int(0) 
maxMoveY = int(0) 
minMoveX = int(0) 
minMoveY = int(0) 
moveBergCount = 0
stuckBergCount = 0
bumpedBergCount = 0
if(advectBergs):
    for i in range(nx-1,0,-1): #i goes from nx-1 to 1 in reverse order
        for j in range(ny-1,0,-1):#j goes backwards as well, as bergs tend to move north
            if(bathymetry[j,i] == 0):
                continue #This is a wall, no bergs here, skip to next iteration
            ## Total freshwater flux / area of bergs in cell 
            ## perhaps this below line could use the actual lambda for the location, but be careful about plume cells
            spd_x = dataMITgcmOcean[2,:,j,i] #[m/s]
            spd_y = dataMITgcmOcean[4,:,j,i] #[m/s]
            for k in range(bergsPerCell[j,i]):
                if(min(bergWidths[k,j,i],bergLength[k,j,i],bergDepths[k,j,i]) > 2):
                    ## Depth check
                    if(bergDepths[k,j,i] > hardMaxDepth):
                        bergDepths[k,j,i] = hardMaxDepth
                        print(f"WARNING deep iceberg adjusted to fit in domain k/j/i: {k}/{j}/{i}")
                    ## assign berg location, advect, re-bin
                    depthIndex = np.max([np.argmin(np.abs(z[:] - bergDepths[k,j,i])),1]) # shallow bergs are entirely in top layer
                    # print(f"depth, ind of berg: {bergDepths[k,j,i]}, {depthIndex}")
                    avgMelt = np.nanmean(dataMITgcm[2,:depthIndex,j,i,],axis=0) #m/day
                    init_x = random.random()*deltaX
                    init_y = random.random()*deltaY
                    ## movement = speed * couplingTimeStep
                    move_x = np.mean(spd_x[:depthIndex]) * couplingTimeStep #only avg to depth of berg
                    move_y = np.mean(spd_y[:depthIndex]) * couplingTimeStep #only avg to depth of berg
                    # print(f"spdx,y {spd_x[:depthIndex]}, {spd_y[:depthIndex]}")
                    advect_x = int(np.floor((init_x + move_x)/deltaX))
                    advect_y = int(np.floor((init_y + move_y)/deltaY))
                    # if(advect_y  < 0):
                    #     print(f"i,j init x,y {i},{j} {init_x:0.2f},{init_y:0.2f}")
                    #     print(f"i,j move x,y {i},{j} {move_x:0.2f},{move_y:0.2f}")
                    #     print(f"i,j advect x,y {i},{j} {advect_x},{advect_y}")
                    if(advect_x > maxMoveX):
                            maxMoveX = advect_x
                    if(advect_y > maxMoveY):
                            maxMoveY = advect_y
                    if(advect_x < minMoveX):
                            minMoveX = advect_x
                    if(advect_y < minMoveY):
                            minMoveY = advect_y
                    if(advect_y != 0 or advect_x != 0):
                        moveBergCount += 1
                    # print(f" j,i {j:2d},{i:2d} advect y, x: {advect_y},{advect_x}")
                    ## move into new cell, increment number of berg in new cell
                    if(i + advect_x < icebergRightHandGate):
                        bergArea = (bergWidths[k,j,i]) * (bergLength[k,j,i])
                        new_j = j + advect_y
                        new_i = i + advect_x
                        ## Bounds on array lengths
                        if(new_j < 0):
                            new_j = 0
                        elif(new_j > ny-1):
                            new_j = ny - 1
                        if(new_i < 0):
                            new_i = 0
                        elif(new_i > nx-1):
                            new_i = nx - 1
                        ## Bounds on where the ocean actually is
                        while(bathymetry[new_j,new_i] == 0):
                            if(new_i < 2): #if inside glacier, push out
                                new_i = 2
                            elif(new_j < ny/2): #if in fjord walls, push back to center
                                new_j = new_j + 1
                            else:
                                new_j = new_j - 1
                        existingAreaHome = np.sum(bergWidthsNew[:,j,i] * bergLengthNew[:,j,i])  #area of bergs in new timestep where you are coming from
                        if((bergArea + existingAreaHome)/(deltaX*deltaY) > maxLambda and new_j == j and new_i == i):
                            ## This means your home cell is full before we even start, and you weren't planning to move
                            ## We are going to say you are now "shoved" by the incoming bergs. You get pushed 1 cell down fjord no matter what. 
                            new_i = i + 1 #bumping you down fjord
                            bumpedBergCount += 1 
                            # print(f"\t\tBumping berg j,i {j},{i} depth {bergDepths[k,j,i]:0.2f}")
                        existingArea = np.sum(bergWidthsNew[:,new_j,new_i] * bergLengthNew[:,new_j,new_i]) #area of bergs in new timestep where you want to go
                        if((bergArea + existingArea)/(deltaX*deltaY) < maxLambda and bergsPerCellNew[new_j,new_i] < maxBergs - 2):
                            bergWidthsNew[bergsPerCellNew[new_j,new_i],new_j,new_i] = bergWidths[k,j,i] - avgMelt * couplingTimeStep/86400
                            bergLengthNew[bergsPerCellNew[new_j,new_i],new_j,new_i] = bergLength[k,j,i] - avgMelt * couplingTimeStep/86400
                            if((bergDepths[k,j,i] - avgMelt * couplingTimeStep/86400) > minBergDepth):
                                bergDepthsNew[bergsPerCellNew[new_j,new_i],new_j,new_i] = bergDepths[k,j,i] - avgMelt * couplingTimeStep/86400
                            else:
                                bergDepthsNew[bergsPerCellNew[new_j,new_i],new_j,new_i] = bergDepths[k,j,i] # too little to melt
                            bergsPerCellNew[new_j,new_i] += 1
                        elif(bergsPerCellNew[j,i] > maxBergs - 2):
                            print(f'**WARNING**  Cell j,i ({j},{i}) is overfull on COUNT, deleting berg depth {bergDepths[k,j,i]:0.2f} m. This is bad.')
                        elif((bergArea + existingAreaHome)/(deltaX*deltaY) > maxLambda):
                            print(f'**WARNING**  Cell j,i ({j},{i}) is overfull on LAMBDA, deleting berg depth {bergDepths[k,j,i]:0.2f} m. This is bad.')
                        else:
                            # print('Cell too full (volume or count), iceberg melts, but does not move')
                            bergWidthsNew[bergsPerCellNew[j,i],j,i] = bergWidths[k,j,i] - avgMelt * couplingTimeStep/86400
                            bergLengthNew[bergsPerCellNew[j,i],j,i] = bergLength[k,j,i] - avgMelt * couplingTimeStep/86400
                            if((bergDepths[k,j,i] - avgMelt * couplingTimeStep/86400) > minBergDepth):
                                bergDepthsNew[bergsPerCellNew[j,i],j,i] = bergDepths[k,j,i] - avgMelt * couplingTimeStep/86400
                            else:
                                bergDepthsNew[bergsPerCellNew[j,i],j,i] = bergDepths[k,j,i] # youre too little!
                            bergsPerCellNew[j,i] += 1
                            stuckBergCount += 1
                    else:
                        # print('Iceberg has left the zone and is lost. Index: %i' %(i + advect))
                        pass
                else: 
                    ## I try to set the minimum berg size such that this doesnt happen in the melange
                    ## You may want to rely on this option like this for freely drifting bergs, but for now its a warning as it should happen
                    print('\t\t** WARNING small iceberg discarded second filter [W,L,D]: [%.3f,%.3f,%.3f] WARNING **' %(bergWidths[k,j,i],bergLength[k,j,i],bergDepths[k,j,i]))
    del i,j,k

    ## Now sorted, we can put the new arrays into the old arrays
    bergWidths = bergWidthsNew.copy()
    bergLength = bergLengthNew.copy()
    bergDepths = bergDepthsNew.copy()
    bergsPerCell = bergsPerCellNew.copy()

print(f'\t\tMax/Min move is (y,x): ({maxMoveY}, {maxMoveX})/({minMoveY}, {minMoveX})')
print('\t\tTotal bergs after move: %i' %np.sum(bergsPerCell))
print(f'\t\tMoving Attempt bergs: {moveBergCount} ({moveBergCount * 100.0 / np.sum(bergsPerCell):0.2f}%)')
print(f'\t\tStuck bergs: {stuckBergCount} ({stuckBergCount * 100.0 / np.sum(bergsPerCell):0.2f}%)')
print(f'\t\tBumped bergs: {bumpedBergCount} ({bumpedBergCount * 100.0 / np.sum(bergsPerCell):0.2f}%)')
if(clipBergs):
    ## Remove Bergs too far along fjord, beyond right hand gate
    bergWidths[:,:,icebergRightHandGate:] = 0
    bergLength[:,:,icebergRightHandGate:] = 0
    bergDepths[:,:,icebergRightHandGate:] = 0
    bergsPerCell[:,icebergRightHandGate:] = 0
print('\t\tTotal bergs after clipping: %i' %np.sum(bergsPerCell))


maxNewBergDepth = 200
maxBergWidth = 0.0642449*maxNewBergDepth**(5/3)
minBergWidth = 40
# We use 2*alpha b/c P(A) ~ A^(-1.9), but A ~ W^2, so P(A) ~ W^(-1.9*2)
alpha = (-1.9)*2
aspectRatio = 1.5

## Making new bergs to fill in the gaps/ ensure left bound stays packing to desired fraction
## importantly we must use length*widths for area as we haven't recalculated OpenFrac
generatedDepths = []
generatedWidths = []
if(addBergs):
    for j in range(ny):
        for i in range(icebergRefreshingGate):
            i += 1 ## we skip the first cell in MITgcm which is the glacier
            if(bathymetry[j,i] == 0):
                continue ## This is a wall, no bergs here
            tmpBergFac = np.sum(bergWidths[:,j,i] * bergLength[:,j,i])/(deltaX*deltaY)
            while(tmpBergFac < icebergCoverLambda):
                ## We keep trying to add while cell is underfull
                randScale = np.random.normal(6,1.22,1)
                randPower = np.random.normal(0.3,0.016,1)
                a = (minBergWidth**(alpha+1) + (maxBergWidth**(alpha+1) - minBergWidth**(alpha+1))*random.random())**(1/(alpha+1))
                b = a / aspectRatio
                if(a > deltaX):
                    a = deltaX * .95 # limit to grid size
                if(b > deltaX):
                    b = deltaX * .95 # Limit to grid size
                if(tmpBergFac + (a*b)/(deltaX*deltaY) < maxLambda): #if this berg is too big, don't get added break, cell left underful
                    bergWidths[bergsPerCell[j,i],j,i] = a
                    bergLength[bergsPerCell[j,i],j,i] = b
                    bergDepths[bergsPerCell[j,i],j,i] = randScale[0] *(a*b)** randPower[0] * (iceDensity/oceanDensity)
                    bergsPerCell[j,i] += 1
                    generatedDepths.append(randScale *(a*b)** randPower * (iceDensity/oceanDensity))
                    generatedWidths.append(a)
                    tmpBergFac = np.sum(bergWidths[:,j,i] * bergLength[:,j,i])/(deltaX*deltaY)
                    # print(np.sum(bergWidthsNew[:,j,i] * bergLengthNew[:,j,i]))
                    # print((deltaX*deltaY))
                    # print(bergsPerCellNew[j,i],"Lambda:%.4f" %tmpBergFac)
                else:
                    # print('Last iceberg was too big, so we exit')
                    break
            # print(tmpBergFac)
    if(len(generatedWidths) == 0):
        generatedDepths.append(0)
        generatedWidths.append(0)

    print('\t\tMax depth/width of new bergs: %.2f/%.2f' %(np.max(generatedDepths),np.max(generatedWidths)))
    plt.figure(3)
    plt.subplot(211)
    plt.hist(np.array(generatedDepths),100)
    plt.xlabel('Depths')
    plt.subplot(212)
    plt.hist(generatedWidths,100)
    plt.xlabel('Widths')
    plt.suptitle('New Berg Distributions')
    plt.tight_layout()
print('\t\tTotal bergs after topping off: %i' %np.sum(bergsPerCell))
# print(np.sum(bergsPerCell,axis=0))


## We now create the new geometry and mask files for MITgcm

## clear these arrays for new version
openFrac[:,:,:] = 1
totalBergArea[:,:,:] = np.nan

for i in range(nx): #loops over whole domain, this calculates new hfill factors
    for j in range(ny):
        numberOfBergs = len(bergDepths[bergDepths[:,j,i] > 0,j,i])
        # print(i,j,numberOfBergs)
        if(numberOfBergs > 0):
            lengths = bergLength[bergDepths[:,j,i] > 0,j,i] #return only non-zeros 
            widths  = bergWidths[bergDepths[:,j,i] > 0,j,i] #return only non-zeros
            depths  = bergDepths[bergDepths[:,j,i] > 0,j,i] #return only non-zeros (negatives should have been blocked before here)
            for k in range(nz):
                cellVolume = deltaX*deltaY*dz[k]
                d_bot = sum_z[k] #bottom of depth bin
                d_top = sum_z[k] - dz[k]
                volume1 = dz[k] * lengths[depths > d_bot] * widths[depths > d_bot]
                SA1 = dz[k]*2*(lengths[depths > d_bot] + widths[depths > d_bot])
                partialFill = (depths < d_bot) & (depths > d_top)
                #partial fill
                volume2 = (depths[partialFill] - d_top) * lengths[partialFill] * widths[partialFill]
                #partial sides
                SA2 = (depths[partialFill] - d_top)*2*(lengths[partialFill] + widths[partialFill]) 
                #bottom
                SA3 = lengths[partialFill] * widths[partialFill]
                #print(np.sum(volume1), np.sum(volume2))
                # print(1-((np.sum(volume1) + np.sum(volume2))/cellVolume))
                openFrac[k,j,i] = 1-((np.sum(volume1) + np.sum(volume2))/cellVolume)
                totalBergArea[k,j,i] = np.sum(SA1) + np.sum(SA2) + np.sum(SA3)
            bergMask[j,i] = 1
            meltMask[j,i] = 1
            ## all free driftin bergs here
            barrierMask[j,i] = 0
            driftMask[j,i] = 1
            if(bergsPerCell[j,i] != numberOfBergs):
                bergsPerCell[j,i] = numberOfBergs
                print('**WARNING** count mis-match of icebergs. Likely you have mis-matched iceberg binary files')
        elif(numberOfBergs == 0):
            openFrac[:,j,i] = 1
            totalBergArea[:,j,i] = 0
            bergMask[j,i] = 0
            meltMask[j,i] = 0
            barrierMask[j,i] = 0
            driftMask[j,i] = 0
del i,j

## Remove blocking from plume locations
barrierMask[plumeMask[:,1] > 1,1] = 0

## Visual outputs, can toggle on/off or save

plt.figure(1)
plt.subplot(221)
varphi = 1-openFrac
varphi[varphi == 1] = np.nan
plt.plot(x[0,:],np.sum(bergsPerCell[1:-1,:],axis=0))
plt.title('berg count')
plt.ylabel('Count [ ]')
plt.xlabel('Along fjord [m]')

plt.subplot(222)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    pc = plt.plot(x[0,:],np.nanmean(varphi[0,:,:],axis=0),linestyle='-')
    pc = plt.plot(x[0,:],np.nanmax(varphi[0,:,:],axis=0),linestyle=':')
# plt.suptitle('$\\lambda$')
plt.ylim([0,1])
plt.ylabel('$\\lambda$')
plt.xlabel('Along fjord [m]')

depthHelper = bergDepths.copy()
depthHelper[depthHelper == 0] = np.nan
plt.subplot(223)
varphi = 1-openFrac
varphi[varphi == 1] = np.nan
# print(varphi[:,:,0:10])
with warnings.catch_warnings():
	warnings.simplefilter("ignore", category=RuntimeWarning)
	pc = plt.pcolormesh(x[0,:],-z,np.nanmean(varphi[:,1:-1,:],axis=1),cmap='gray_r')
cbar = plt.colorbar(pc)
depthHelper = bergDepths.copy()
depthHelper[depthHelper == 0] = np.nan
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    plt.plot(x[0,:],-np.nanmean(depthHelper[:,1:-1,:],axis=(0,1)),color='xkcd:orange',linestyle='--',label='$<D>$')
plt.suptitle('$\\varphi$')
plt.legend()
plt.ylabel('Depth [m]')
plt.xlabel('Along fjord [m]')
cbar.set_label('ice fraction')

plt.subplot(224)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    pc = plt.pcolormesh(x[0,:],y[:,0],varphi[0,:,:],cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
plt.suptitle('$\\lambda$')
plt.ylabel('Across fjord [m]')
plt.xlabel('Along fjord [m]')
cbar.set_label('ice fraction')
plt.clim([0,1])
del depthHelper

print('\t\tMax lambda/Depth/Width is: %.3f/%.3f/%.3f' %(np.max(varphi[0,:,:]),np.max(bergDepths),np.max(bergWidths)))

plt.tight_layout() 
## you can save this figure if you'd like, it is a fun movie if you combine them all
plt.savefig('figs/advectBergs%08d.png' %maxStep,format='png',dpi=150)

## I toggle this for troubleshooting purposes
if(not trialOnly):
	# Write files that need updating
	write_bin('bergMask.bin',bergMask)
	# write_bin('bergMaskNums.bin',bergMaskNums)
	write_bin('numBergsPerCell.bin',bergsPerCell)
	write_bin('openFrac.bin',openFrac)
	write_bin('totalBergArea.bin',totalBergArea)
	write_bin('meltMask.bin',meltMask)
	write_bin('driftMask.bin',driftMask)
	write_bin('barrierMask.bin',barrierMask)
	write_bin('icebergs_depths.bin',bergDepths)
	write_bin('icebergs_widths.bin',bergWidths)
	write_bin('icebergs_length.bin',bergLength)
else:
    print('Not Saving')
    print('Seconds to advect: %.4f' % (time.time() - start_time))
if(showPlots):
    plt.show()
    plt.close()
# os.rename("totalBergArea.bin", "totalBergArea.bin" + "~") #save old file incase something goes wrong
# os.rename("openFrac.bin", "openFrac.bin" + "~") #save old file incase something goes wrong
