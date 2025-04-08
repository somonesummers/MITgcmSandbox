import os
import sys
import numpy as np
import fileinput
from MITgcmutils import mds
import matplotlib.pyplot as plt
import warnings
from scipy.interpolate import interp1d
import cmocean
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import constants, glaciome
import glob
import pickle
import time
import random

start_time = time.time()

def write_bin(fname, data):
    # print(fname + " " + str(np.shape(data)))
    data.astype(">f8").tofile('input/'+fname)

# Controllers for which steps of advection process we do

advectBergs = True
clipBergs = True
shadedEdge = True
fillIn = True
# icebergCoverLambda = .80 
couplingTimeStep = 1*24*3600 #[s]

iceDensity = 917 # [kg/m^3]
oceanDensity = 1030 # [kg/m^3]

# Some rules on advecting/melting/sizing
minBergDepth = 5 # bergs not allowed to shrink/melt smaller than this [m]
maxLambda = .80 # Cells are full at this lambda
scaleMin = 0.2 # percent error allowed between MITgcm effective depth and GLACIOME depth
# Bergs added to top off to icebergCoverLamba left of refreshGate
icebergRefreshingGate = 3
maxBergs = 500


# Load model wide parameters

sys.path.append('.')
from melangeModel import *

# Get most recent data from MITgcm

maxStep = 0
for file in os.listdir('results'):
    # print(file)
    if "BRGFlx.0" in file:
        words = file.split(".")
        # print(words[1])  
        if int(words[1]) > maxStep:
            maxStep = int(words[1])

# Iceberg melt data
dataMITgcm = mds.rdmds("results/BRGFlx", maxStep)

# Geometry of MITgcm domain
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

#This needs to know about z. Bergs cant be deeper than this
hardMaxDepth = np.max(np.abs(z)) - 5
#icebergs are destroyed when they move past here. Issues when GLACIOME grows past here
icebergRightHandGate = int(nx*0.9)

# get most recent Glaciome1D data. This drives advection speeds and geometry updates
files = sorted(glob.glob('couplingResults/MITgcmRun_[0-9][0-9][0-9][0-9][0-9].pickle'))
# print(files)
with open(files[-1], 'rb') as file:
    data = pickle.load(file)
    file.close()
    index = int(str(file).split(".")[1].split("_")[1])

# We load interpolators to apply to MITgcm grid
U_glaciome1d = (data.U+data.Ut-data.Uc)/(24*3600*365) #raw is m/yr, convert to m/s
H = np.concatenate(([data.H0], data.H, [data.HL]))
X_glaciome1d = data.X + deltaX #meters, shift to align with front of MITgcm melange
X_ = np.concatenate(([data.X[0]],data.X_,[data.X[-1]])) + deltaX
glaciome_L = data.L + deltaX
speedLookup = interp1d(X_glaciome1d,U_glaciome1d,bounds_error=False,fill_value=[0])
thicknessLookup = interp1d(X_,H,bounds_error=False,fill_value=[0])
# plt.figure
# plt.plot(X_glaciome1d,U_glaciome1d)
# plt.show()

# Now we load the iceberg geometry and mask files

#2D Fields 
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
#3D Fields
totalBergArea = np.fromfile('input/totalBergArea.bin', dtype='>f8')
totalBergArea = totalBergArea.reshape((nz,ny,nx))
openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
openFrac = openFrac.reshape((nz,ny,nx))

# lam = np.nanmean((1-openFrac[0,1:-1,:]),axis=0)
# lam2 = np.mean(np.nansum(bergDepths[:,1:-1,:]*bergWidths[:,1:-1,:]*bergLength[:,1:-1,:],axis=0)/(deltaX*deltaY),axis=0)
# print('lambda',lam,lam2)

# Visual midpoints, can toggle on/off or save

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
plt.ylim([0,0.9])
plt.ylabel('$\\lambda$')
plt.xlabel('Along fjord [m]')

# Create fresh of geometry arrays, this lets us move bergs without double counting

bergWidthsNew   = np.zeros_like(bergWidths)
bergLengthNew   = np.zeros_like(bergWidths)
bergDepthsNew   = np.zeros_like(bergWidths)
bergsPerCellNew = np.zeros_like(bergsPerCell).astype(np.int64)

print('Total bergs before move: %i' %np.sum(bergsPerCell))

## Advect bergs. We go backwards (front to back) to ensure overfull icebergs can get priority to stay where they are
# If you want it to not really be random
# np.random.seed(3)
maxMove = int(0) #track this to know how far the farthest berg moved. We only 'top off' this far.
if(advectBergs):
    for j in range(ny):
        for i in range(nx-1,0,-1): #i goes from nx-1 to 1 in reverse order
                # We can skip the first cell as that is a glacier in MITgcm

                # Total freshwater flux / area of bergs in cell 
                # perhaps this below line could use the actual lambda for the location, but be careful about plume cells
                effectiveMelt = np.nanmean(np.nansum(dataMITgcm[0,:,:,1:],axis=0),axis=0)/(deltaX*deltaY*icebergCoverLambda)*(24*3600)*1000/917 #[m]
                for k in range(bergsPerCell[j,i]):
                    if(min(bergWidths[k,j,i],bergLength[k,j,i],bergDepths[k,j,i]) > 5):
                        #assign berg location, advect, re-bin
                        init_x = random.random()*deltaX
                        movement = speedLookup(x[0,i]) * couplingTimeStep
                        # movement = speed * couplingTimeStep
                        # print(movement)
                        advect = int(round((init_x + movement)/deltaX))
                        if(i == 1 and advect > maxMove):
                                maxMove = advect
                        depthIndex = np.argmin(np.abs(z[:]- bergDepths[k,j,i]))
                        # print(advect)
                        #move into new cell, increment number of berg in new cell
                        if(i + advect < icebergRightHandGate):
                            bergArea = (bergWidths[k,j,i]) * (bergLength[k,j,i])
                            existingArea = np.sum(bergWidthsNew[:,j,i + advect] * bergLengthNew[:,j,i + advect])
                            if((bergArea + existingArea)/(deltaX*deltaY) < maxLambda):
                                # print('Moved,',(bergArea + existingArea)/(deltaX*deltaY))
                                bergWidthsNew[bergsPerCellNew[j,i + advect],j,i + advect] = bergWidths[k,j,i] #- avgMelt[depthIndex]
                                bergLengthNew[bergsPerCellNew[j,i + advect],j,i + advect] = bergLength[k,j,i] #- avgMelt[depthIndex]
                                if((bergDepths[k,j,i] - effectiveMelt[i]) > minBergDepth):
                                    bergDepthsNew[bergsPerCellNew[j,i + advect],j,i + advect] = bergDepths[k,j,i] - effectiveMelt[i]
                                else:
                                    bergDepthsNew[bergsPerCellNew[j,i + advect],j,i + advect] = bergDepths[k,j,i]
                                # if(avgMelt[depthIndex] > 0):
                                    # print('%f %f %f to %f %f %f' %(bergWidths[k,j,i],bergLength[k,j,i],bergDepths[k,j,i],
                                    #     bergWidths[k,j,i]-avgMelt[depthIndex],
                                    #     bergLength[k,j,i]-avgMelt[depthIndex],
                                    #     bergDepths[k,j,i]-avgMelt[depthIndex]))
                                bergsPerCellNew[j,i + advect] += 1
                            else:
                                # print('Cell too full, iceberg melts, but does not move')
                                bergWidthsNew[bergsPerCellNew[j,i],j,i] = bergWidths[k,j,i]# - avgMelt[depthIndex]
                                bergLengthNew[bergsPerCellNew[j,i],j,i] = bergLength[k,j,i]# - avgMelt[depthIndex]
                                bergDepthsNew[bergsPerCellNew[j,i],j,i] = bergDepths[k,j,i] - effectiveMelt[i]
                                bergsPerCellNew[j,i] += 1
                        else:
                            # print('Iceberg has left the zone and is lost. Index: %i' %(i + advect))
                            pass
                    else:
                        print('small iceberg discarded second filter [W,L,D]:', ([bergWidths[k,j,i],bergLength[k,j,i],bergDepths[k,j,i]]))
    del i,j,k

    # Now sorted, we can put the new arrays into the old arrays
    bergWidths = bergWidthsNew.copy()
    bergLength = bergLengthNew.copy()
    bergDepths = bergDepthsNew.copy()
    bergsPerCell = bergsPerCellNew.copy()

#Track where the edge of the melange is according to GLACIOME
melangeIndex = int(np.ceil(glaciome_L/deltaX))
print('Max move is',maxMove)
print('Total bergs after move: %i' %np.sum(bergsPerCell))

if(clipBergs):
    # Remove Bergs too far along fjord
    bergWidths[:,:,melangeIndex:] = 0
    bergLength[:,:,melangeIndex:] = 0
    bergDepths[:,:,melangeIndex:] = 0
    bergsPerCell[:,melangeIndex:] = 0
print('Total bergs after clipping: %i' %np.sum(bergsPerCell))

if(shadedEdge):
    # Trim final cell, this is the last cell melange extends into
    # Should trim proportional to how far into this cell the melange extends. Clipping by number of bergs remains somewhat random
    trimFac = (glaciome_L - (melangeIndex-1)*deltaX)/deltaX
    bergsPerCell[:,melangeIndex-1] = (bergsPerCell[:,melangeIndex-1]*trimFac).astype(np.int64)
    for j in range(ny):
        # print('j = %i, k = %i' %(j,bergsPerCellNew[j,melangeIndex-1]))
        bergWidths[bergsPerCell[j,melangeIndex-1]:,j,melangeIndex-1] = 0
        bergLength[bergsPerCell[j,melangeIndex-1]:,j,melangeIndex-1] = 0
        bergDepths[bergsPerCell[j,melangeIndex-1]:,j,melangeIndex-1] = 0
print('Total bergs after shadedEdge: %i' %np.sum(bergsPerCell))

if(fillIn):
    # We fill in the backside of the mélange with new icebergs breaking off the glacier.
    # print(np.sum(bergsPerCellNew,axis=0))
    maxBergDepth = np.min([data.H0*2.5, z[-1]*.8])
    maxBergWidth = 0.0642449*maxBergDepth**(5/3)
    minBergWidth = 40
    # We use 2*alpha b/c P(A) ~ A^(-1.9), but A ~ W^2, so P(A) ~ W^(-1.9*2)
    alpha = (-1.9)*2
    aspectRatio = 1.5

    ## Making new bergs to fill in the gaps/ ensure left bound stays packing to desired fraction
    # importantly we must use length*widths for area as we haven't recalculated OpenFrac
    print('Max depth/width of new bergs: %.2f, %.2f' %(maxBergDepth,maxBergWidth))
    generatedDepths = []
    generatedWidths = []
    for j in range(ny-2):
        j += 1 # we skip the fjord walls
        for i in range(np.max([maxMove,icebergRefreshingGate])):
            i += 1 # we skip the first cell in MITgcm which is the glacier
            tmpBergFac = np.sum(bergWidths[:,j,i] * bergLength[:,j,i])/(deltaX*deltaY)
            while(tmpBergFac < icebergCoverLambda):
                # We keep trying to add while cell is underfull
                randScale = np.random.normal(6,1.22,1)
                randPower = np.random.normal(0.3,0.016,1)
                a = (minBergWidth**(alpha+1) + (maxBergWidth**(alpha+1) - minBergWidth**(alpha+1))*random.random())**(1/(alpha+1))
                b = a / aspectRatio
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

    plt.figure(3)
    plt.subplot(211)
    plt.hist(np.array(generatedDepths),100)
    plt.xlabel('Depths')
    plt.subplot(212)
    plt.hist(generatedWidths,100)
    plt.xlabel('Widths')
    print('Total bergs after topping off: %i' %np.sum(bergsPerCell))
# print(np.sum(bergsPerCell,axis=0))

## So now we have all our old icebergs added, edge shaded, new bergs added.
# We now create the new geometry and mask files for MITgcm

#nuke these arrays for new version
openFrac[:,:,:] = 1
totalBergArea[:,:,:] = np.nan

# We need to rescale our bergs according to the GLACIOME1D geometry
# We get effective depth from the geometry files (how deep would all this ice be if packed 100%)
effectiveDepth = np.nanmean(np.nansum(bergDepths[:,1:-1,:]*bergWidths[:,1:-1,:]*bergLength[:,1:-1,:],axis=0)/(deltaX*deltaY),axis=0)

print('Melange edge at i=%i' %melangeIndex)
for i in range(nx): #don't scale up the shaded edge
    for j in range(ny):
        numberOfBergs = len(bergDepths[bergDepths[:,j,i] > 0,j,i])
        # print(i,j,numberOfBergs)
        if(numberOfBergs > 0):
            if(i < (melangeIndex - 1)): #Only scale for bergs within melange, not in shadedEdge
                scaleFactor =  thicknessLookup(x[0,i]) * iceDensity/oceanDensity * icebergCoverLambda / effectiveDepth[i]
                scaleFactor = (scaleFactor if abs(scaleFactor - 1) > scaleMin else 1)
                if(scaleFactor != 1 and j == 1):
                    print('bergs are adjusted: %.3f i=%i, nB: %i' %(scaleFactor,i,np.sum(bergsPerCell[:,i])))
                bergDepths[:,j,i] = bergDepths[:,j,i] * scaleFactor
                bergDepths[bergDepths[:,j,i] > hardMaxDepth,j,i] = hardMaxDepth
                bergDepths[(bergDepths[:,j,i] > 0) & (bergDepths[:,j,i] < minBergDepth),j,i] = minBergDepth + .01 # set min depth, bump 'em up a bit
            lengths = bergLength[bergLength[:,j,i] > 0,j,i] #return only non-zeros
            widths = bergWidths[bergWidths[:,j,i] > 0,j,i] #return only non-zeros
            depths = bergDepths[bergDepths[:,j,i] > 0,j,i] #return only non-zeros
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
            barrierMask[j,i] = 1
            driftMask[j,i] = 0
            if(bergsPerCell[j,i] != numberOfBergs):
                bergsPerCell[j,i] = numberOfBergs
                print('**WARNING** count mis-match')
        elif(numberOfBergs == 0):
            openFrac[:,j,i] = 1
            totalBergArea[:,j,i] = 0
            bergMask[j,i] = 0
            meltMask[j,i] = 0
            barrierMask[j,i] = 0
            driftMask[j,i] = 0
del i,j

#remove blocking from plume locations
barrierMask[plumeMask[:,1] > 1,1] = 0

plt.figure(1)
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
    plt.plot(x[0,:],-effectiveDepth/icebergCoverLambda,color='xkcd:red',linestyle='--')
    plt.plot(x[0,:],-np.nanmean(depthHelper[:,1:-1,:],axis=(0,1)),color='xkcd:orange',linestyle='--')
plt.plot(X_,-H*iceDensity/oceanDensity,color='xkcd:green',linestyle='--')
plt.plot([glaciome_L,glaciome_L],[20,-600],linestyle='--',color='cyan')
plt.suptitle('$\\varphi$')
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

print('Max lambda is:',np.max(varphi[0,:,:]))

plt.tight_layout() 
plt.savefig('figs/advectBergs%08d.png' %maxStep,format='png',dpi=150)

# I toggle this for troubleshooting purposes
if(True):
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
    plt.show()
    plt.close()
# os.rename("totalBergArea.bin", "totalBergArea.bin" + "~") #save old file incase something goes wrong
# os.rename("openFrac.bin", "openFrac.bin" + "~") #save old file incase something goes wrong