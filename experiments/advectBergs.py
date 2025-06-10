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
icebergFringe = False
trialAndPlot = False
couplingTimeStep = 1*24*3600 #[s]

iceDensity = 917 # [kg/m^3]
oceanDensity = 1030 # [kg/m^3]

# Some rules on advecting/melting/sizing
minBergDepth = 5 # bergs not allowed to shrink/melt smaller than this [m]
maxBergDepth = 450 # bergs not allowed to shrink/melt smaller than this [m]
maxLambda = .90 # Cells are full at this lambda [hard limit]
scaleMin = 0.2 # percent error allowed between MITgcm effective depth and GLACIOME depth
# Bergs added to top off to icebergCoverLamba left of refreshGate
icebergRefreshingGate = 3
maxBergs = 500


# Load model wide parameters

sys.path.append('.')
from melangeModel import *

icebergPhi = np.min([.8,icebergCoverLambda]) # Scaling from GLACIOME1D
icebergCoverLambda = np.min([.9,icebergCoverLambda * 2]) # How full filled cells are

print('\t\tLambda/Phi %.3f/%.3f' %(icebergCoverLambda,icebergPhi))
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
#3D Fields
totalBergArea = np.fromfile('input/totalBergArea.bin', dtype='>f8')
totalBergArea = totalBergArea.reshape((nz,ny,nx))
openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
openFrac = openFrac.reshape((nz,ny,nx))

# lam = np.nanmean((1-openFrac[0,1:-1,:]),axis=0)
# lam2 = np.mean(np.nansum(bergDepths[:,1:-1,:]*bergWidths[:,1:-1,:]*bergLength[:,1:-1,:],axis=0)/(deltaX*deltaY),axis=0)
# print('lambda',lam,lam2)

# Create fresh of geometry arrays, this lets us move bergs without double counting

bergWidthsNew   = np.zeros_like(bergWidths)
bergLengthNew   = np.zeros_like(bergWidths)
bergDepthsNew   = np.zeros_like(bergWidths)
bergsPerCellNew = np.zeros_like(bergsPerCell).astype(np.int64)

print('\t\tTotal bergs at start move: %i' %np.sum(bergsPerCell))

#Track where the edge of the melange is according to GLACIOME
melangeIndex = int(np.ceil(glaciome_L/deltaX))
print(f'\t\tMélange edge is index: {melangeIndex}')
if(icebergFringe): 
    # Remove Bergs beyong the melange, ie the fringe. 
    # This must bedone before advecting mélange bergs to ensure no big icebergs get 'stuck' in the mélange.
    # Once iceberg fringe advection is enabled, this may be possible to remove.
    bergWidths[:,:,melangeIndex:] = 0
    bergLength[:,:,melangeIndex:] = 0
    bergDepths[:,:,melangeIndex:] = 0
    bergsPerCell[:,melangeIndex:] = 0
print('\t\tTotal bergs after first clipping, before move: %i' %np.sum(bergsPerCell))

## Advect bergs. We go backwards (front to back) to ensure overfull icebergs can get priority to stay where they are
# If you want it to not really be random
# np.random.seed(3)
maxMove = int(0) #track this to know how far the farthest berg moved. We only 'top off' this far.
if(advectBergs):
    for j in range(ny):
        for i in range(nx-1,0,-1): #i goes from nx-1 to 1 in reverse order
                if(bathymetry[j,i] == 0):
                    continue #This is a wall, no bergs here
                # Total freshwater flux / area of bergs in cell 
                # perhaps this below line could use the actual lambda for the location, but be careful about plume cells
                effectiveMelt = np.nanmean(np.nansum(dataMITgcm[0,:,:,1:],axis=0),axis=0)/(deltaX*deltaY*icebergPhi)*(24*3600)*1000/917 #[m]
                for k in range(bergsPerCell[j,i]):
                    if(min(bergWidths[k,j,i],bergLength[k,j,i],bergDepths[k,j,i]) > 2):
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
                                    bergDepthsNew[bergsPerCellNew[j,i + advect],j,i + advect] = bergDepths[k,j,i] - effectiveMelt[i] # melt big bergs
                                else:
                                    bergDepthsNew[bergsPerCellNew[j,i + advect],j,i + advect] = bergDepths[k,j,i] # too little to melt
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
                                if((bergDepths[k,j,i] - effectiveMelt[i]) > minBergDepth):
                                    bergDepthsNew[bergsPerCellNew[j,i],j,i] = bergDepths[k,j,i] - effectiveMelt[i] # melt big bergs
                                else:
                                    bergDepthsNew[bergsPerCellNew[j,i],j,i] = bergDepths[k,j,i] # youre too little!
                                bergsPerCellNew[j,i] += 1
                        else:
                            # print('Iceberg has left the zone and is lost. Index: %i' %(i + advect))
                            pass
                    else: 
                        # I try to set the minimum berg size such that this doesnt happen in the melange
                        # You may want to rely on this option like this for freely drifting bergs, but for now its a warning as it should happen
                        print('\t\t** WARNING small iceberg discarded second filter [W,L,D]: [%.3f,%.3f,%.3f] WARNING **' %(bergWidths[k,j,i],bergLength[k,j,i],bergDepths[k,j,i]))
    del i,j,k

    # Now sorted, we can put the new arrays into the old arrays
    bergWidths = bergWidthsNew.copy()
    bergLength = bergLengthNew.copy()
    bergDepths = bergDepthsNew.copy()
    bergsPerCell = bergsPerCellNew.copy()

print('\t\tMax move is',maxMove)
print('\t\tTotal bergs after move: %i' %np.sum(bergsPerCell))

if(clipBergs):
    # Remove Bergs too far along fjord
    bergWidths[:,:,melangeIndex:] = 0
    bergLength[:,:,melangeIndex:] = 0
    bergDepths[:,:,melangeIndex:] = 0
    bergsPerCell[:,melangeIndex:] = 0
print('\t\tTotal bergs after clipping: %i' %np.sum(bergsPerCell))

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
print('\t\tTotal bergs after shadedEdge: %i' %np.sum(bergsPerCell))

if(fillIn):
    # We fill in the backside of the mélange with new icebergs breaking off the glacier.
    # print(np.sum(bergsPerCellNew,axis=0))
    maxNewBergDepth = np.min([data.H0*2.5, z[-1]*.8])
    maxBergWidth = 0.0642449*maxNewBergDepth**(5/3)
    minBergWidth = 40
    # We use 2*alpha b/c P(A) ~ A^(-1.9), but A ~ W^2, so P(A) ~ W^(-1.9*2)
    alpha = (-1.9)*2
    aspectRatio = 1.5

    ## Making new bergs to fill in the gaps/ ensure left bound stays packing to desired fraction
    # importantly we must use length*widths for area as we haven't recalculated OpenFrac
    generatedDepths = []
    generatedWidths = []
    for j in range(ny):
        for i in range(np.max([maxMove,icebergRefreshingGate])):
            i += 1 # we skip the first cell in MITgcm which is the glacier
            if(bathymetry[j,i] == 0):
                continue #This is a wall, no bergs here
            tmpBergFac = np.sum(bergWidths[:,j,i] * bergLength[:,j,i])/(deltaX*deltaY)
            while(tmpBergFac < icebergCoverLambda):
                # We keep trying to add while cell is underfull
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

## So now we have all our old icebergs added, edge shaded, new bergs added.
# We now create the new geometry and mask files for MITgcm

#nuke these arrays for new version
openFrac[:,:,:] = 1
totalBergArea[:,:,:] = np.nan

# We need to rescale our bergs according to the GLACIOME1D geometry
# We get effective depth from the geometry files (how deep would all this ice be if packed 100%)
effectiveDepth = np.nanmean(np.nansum(bergDepths[:,:,:]*bergWidths[:,:,:]*bergLength[:,:,:],axis=0)/(deltaX*deltaY),axis=0)

if(icebergFringe):
    from bisect import bisect_left
    def find_closest_indices(sorted_A, sorted_B):
        closest_indices = []
        for a in sorted_A:
            pos = bisect_left(sorted_B, a)  # Find position in B where a would fit
            # Compare neighbors to find the closest
            if pos == 0:
                closest_indices.append(0)
            elif pos == len(sorted_B):
                closest_indices.append(len(sorted_B) - 1)
            else:
                before = pos - 1
                after = pos
                closest_indices.append(before if abs(sorted_B[before] - a) <= abs(sorted_B[after] - a) else after)
        return closest_indices
    ## Making new bergs to fill ocean infront of mélange to desired fill (20 - 1 percent over 10km)
    ## We generate these in bulk then assign them to positions
    
    #Geometry of the fringe
    bergType = 1 # 1 = block 2 = cone (not implemented)
    scaling = 1 # 1 = Sulak 2017 2 = Barker 2004
    fringeMaxBergDepth = 160 # (m) 
    fringeLengthInd = int(20000/deltaX)
    print(f'\t\tIceberg fringe is {fringeLengthInd} cells long')
    fringeBergConc = np.zeros(np.shape(bergDepths[0,:,:]))
    fringeBergMask = np.zeros_like(bergMask)
    if(melangeIndex + fringeLengthInd > nx):
        print('Iceberg Fringe extends beyond domain, errors to follow')
    
    #set concentration (lambda) of fringe
    fringeBergConc[:,melangeIndex:(melangeIndex+fringeLengthInd)] = np.linspace(20,5,num=fringeLengthInd)
    fringeBergConc[bathymetry == 0] = 0 #walls get no icebergs
    fringeBergMask[fringeBergConc > 0] = 1

    desiredBergArea = np.sum(fringeBergConc/100.0*deltaX*deltaY)
    bergMaskArea = np.sum(fringeBergMask*deltaX*deltaY)
    print('Area where bergs live: ' + str(bergMaskArea) + ' m^2')
    print('Desired berg area: ' + str(desiredBergArea) + ' m^2')
    print('Ratio: ' + str(desiredBergArea/bergMaskArea*100) + '%')

    if(scaling == 1): # then use Sulak17 volume-area scaling volume = 6.0*area^1.30
        # assumes volume = L*W*D and W = L/1.62 (Dowdeswell et al 1992)
        if(maxBergWidth==0):
            maxBergWidth = 0.0642449*fringeMaxBergDepth**(5/3)
            # minBergWidth = 0.0642449*minBergDepth**(5/3)
        elif(fringeMaxBergDepth==0):
            fringeMaxBergDepth = 5.19155*maxBergWidth**(5/3)
            minBergDepth = 5.19155*minBergWidth**(5/3)
    elif(scaling == 2): # Then use Barker04 width-depth relationship
        # Depth = 2.91*Width^0.71
        if(maxBergWidth==0):
            maxBergWidth = (100*10**(58/71)*fringeMaxBergDepth**(100/71)) / (291*291**(29/71))
            #minBergWidth = (100*10**(58/71)*minBergDepth**(100/71)) / (291*291**(29/71))        
        elif(fringeMaxBergDepth==0):
            fringeMaxBergDepth = 2.91*maxBergWidth^0.71
            minBergDepth = 2.91*minBergWidth^0.71

    numberOfBergs = 25 #low start, immediately doubled by scheme below, so guess low, high guesses (300%+) can cause to fail
    bergTopArea = 0
    areaResidual = 1
    # Generate the Inverse Power Law cumulative distribution function
    # over the range minBergWidth-maxBergWidth with a slope of alpha.
    print('Making bergs, this can take a few loops...')
    loop_count = 1

    # np.random.seed(2)
    # print('random seed set, not really random anymore')

    while(np.abs(areaResidual) > .005 ): # Create random power dist of bergs, ensure correct surface area
        numberOfBergs = round(numberOfBergs * (1 + areaResidual))  
        print('\tnumberOfBergs: ' + str(numberOfBergs))
        x_width = np.arange(minBergWidth, maxBergWidth, (maxBergWidth-minBergWidth)/(numberOfBergs*1e2))
        inversePowerLawPDF_width = ((-alpha-1) / minBergWidth) * (x_width/minBergWidth) ** (alpha)
            # Get the CDF numerically
        inversePowerLawCDF_width = np.cumsum(inversePowerLawPDF_width)
            # Normalize
        inversePowerLawCDF_width = inversePowerLawCDF_width / inversePowerLawCDF_width[-1]
            
            # Generate number_of_bergs uniformly distributed random numbers.
        uniformlyDistributedRandomNumbers = np.random.uniform(0,1,numberOfBergs)
        
        inversePowerLawDistNumbers_width = np.zeros(uniformlyDistributedRandomNumbers.size);

        nearestIndex_width = [0] * uniformlyDistributedRandomNumbers.size

        nearestIndex_width = find_closest_indices(uniformlyDistributedRandomNumbers,inversePowerLawCDF_width)


        inversePowerLawDistNumbers_width = x_width[nearestIndex_width];
        inversePowerLawDistNumbers_length = inversePowerLawDistNumbers_width/1.12 # Widths are bigger 

        randScale = np.random.normal(6,1.22,numberOfBergs)
        randPower = np.random.normal(0.3,0.016,numberOfBergs)
        inversePowerLawDistNumbers_depth = randScale * (inversePowerLawDistNumbers_width*inversePowerLawDistNumbers_length) ** randPower * (920/1025)
        inversePowerLawDistNumbers_depth[inversePowerLawDistNumbers_depth < 5.123] = 5.123 #doest like round numbers, cap low end of bergs

        tooWide = np.count_nonzero(inversePowerLawDistNumbers_width > deltaX * np.sqrt(maxLambda))  #disallow completely full cells
        tooLong = np.count_nonzero(inversePowerLawDistNumbers_length > deltaX * np.sqrt(maxLambda))
        inversePowerLawDistNumbers_width[inversePowerLawDistNumbers_width > deltaX * np.sqrt(maxLambda)] = deltaX * np.sqrt(maxLambda) - .01  # Max width is grid cell (assumed square)
        inversePowerLawDistNumbers_length[inversePowerLawDistNumbers_length > deltaX * np.sqrt(maxLambda)] = deltaX * np.sqrt(maxLambda) - .01 # Max length is grid cell (assumed square)
        if(tooLong + tooWide > 0):
            print('\t\tBergs clipped: %i for width, %i for length' % (tooWide, tooLong))
        
        bergTopArea = sum(inversePowerLawDistNumbers_width*inversePowerLawDistNumbers_length)
        areaResidual = (desiredBergArea - bergTopArea)/desiredBergArea
        print('\t\t%.2f %% Bergs' % (bergTopArea/bergMaskArea*100))
        print('\t\tareaResidual %.2f %%' % (areaResidual * 100))
        loop_count += 1
    print('====== Success! Found our bergs =====')
    print('Width min/mean/max: %f/%f/%f [m]' % (np.min(inversePowerLawDistNumbers_width),np.mean(inversePowerLawDistNumbers_width),np.max(inversePowerLawDistNumbers_width)))
    print('Depth min/mean/max: %f/%f/%f [m]' % (np.min(inversePowerLawDistNumbers_depth),np.mean(inversePowerLawDistNumbers_depth),np.max(inversePowerLawDistNumbers_depth)))
    print('Total Berg Area %f' % bergTopArea)
    print('Total Berg fract: %.2f %%' % (bergTopArea/bergMaskArea*100))

    # Now we sort these berg into cell, randomly 
    bergMaski = 0  #Bad name, but this is the count of cells that will recieve bergs
    bergDict = {}

    for j in range(ny):
        for i in range(nx):
            if(fringeBergMask[j,i] == 1):
                # print('i,j, bergmask',i,j,fringeBergMask[j,i])
                bergMaski = 1 + bergMaski #Needs to start at 1, as non-bergs will be 0
                bergMaskNums[j,i] = bergMaski #Assign Mask Nums, not random as we'll randomly place bergs in cells
                bergDict[bergMaski] = [j,i] #This lets us do 1-D loops for the whole grid
    print('%i cells with bergs' % bergMaski)

    # Sort my bergs
    sorted_indices = np.argsort(-inversePowerLawDistNumbers_depth) # Sort backwards to get descending from big to small bergs 
    sorted_depth = inversePowerLawDistNumbers_depth[sorted_indices]
    sorted_width = inversePowerLawDistNumbers_width[sorted_indices]
    sorted_length = inversePowerLawDistNumbers_length[sorted_indices]
    assignedCell = np.random.randint(0,bergMaski,[numberOfBergs]) # In this script, every berg has a home

    # Array for bergs
    bergsPerCellLimit = 500
    icebergs_depths = np.zeros([bergMaski,bergsPerCellLimit])
    icebergs_widths = np.zeros([bergMaski,bergsPerCellLimit])
    icebergs_length = np.zeros([bergMaski,bergsPerCellLimit])  #careful, not plural as to length match

    # np.random.seed(2)
    assignedCell = np.random.randint(0,bergMaski,[numberOfBergs]) #every Berg has a spot

    icebergs_per_cell = np.zeros([bergMaski],dtype=np.int16)
    icebergs_area_per_cell = np.zeros([bergMaski])

    for i in range(numberOfBergs): 
        j = assignedCell[i]
        # print('looking at mask number',j,'at berg',i)
        # print('Berg number', icebergs_per_cell[j],'in this cell')
        # print('Width, Length',sorted_width[i], sorted_length[i])
        bergArea = sorted_width[i] * sorted_length[i]
        loopLimiter = 0
        while(bergArea > (deltaX * deltaY  * (fringeBergConc[bergDict[j+1][0],bergDict[j+1][1]]/100) - icebergs_area_per_cell[j])): #if above 'full', pick random new cell, accept with decreasing probability
            j_old = j
            j = np.random.randint(0,bergMaski)
            loopLimiter += 1
            if((bergArea + icebergs_area_per_cell[j])/(deltaX * deltaY) < maxLambda - .01): #only consider accepting if under 95
                odds = np.abs(np.random.normal(0,.5,1))  #randomly accepts those that are big in overfull cells, but at decreasing frequency
                overFull = ((bergArea + icebergs_area_per_cell[j])/(deltaX * deltaY)*100 - fringeBergConc[bergDict[j+1][0],bergDict[j+1][1]])
                if(odds > overFull):
                    # print('accepting overfull')
                    assignedCell[i] = j  #if we it a shuffling critera
                    break
            if(loopLimiter > bergMaski*20): #eventually we have to force some in
                indexesAllowed = np.where((deltaX * deltaY * maxLambda - icebergs_area_per_cell)  > bergArea)
                # print(indexesAllowed)
                randi = np.random.randint(0,len(indexesAllowed[0]))
                j = indexesAllowed[0][randi]
                assignedCell[i] = j  #if we it a shuffling critera, must line up for calculation below
                print('\t Randomly missed, will force into cell with room: %i' % j)
                if((np.min(icebergs_area_per_cell) + bergArea)/(deltaX * deltaY) > maxLambda):
                    print('WARNING cell very full: %.2f%%' %((np.min(icebergs_area_per_cell) + bergArea)*100/(deltaX * deltaY)))
                break
        
        icebergs_depths[j,icebergs_per_cell[j]] = sorted_depth[i]
        icebergs_widths[j,icebergs_per_cell[j]] = sorted_width[i]
        icebergs_length[j,icebergs_per_cell[j]] = sorted_length[i]
        icebergs_per_cell[j] += 1
        # icebergs_area_per_cell[j] = np.sum(icebergs_widths[j,:]*icebergs_length[j,:])
        icebergs_area_per_cell[j] += bergArea
    #print('Bergs per cell and filled faction at surface for spot check')     
    #print(icebergs_per_cell)
    #print(np.round(icebergs_area_per_cell/(deltaX*deltaY),2))
    print('Max fill for fringe is: %.2f%%' % (np.nanmax(icebergs_area_per_cell/(deltaX*deltaY))*100))

    #above makes 1D list of all bergs using a dictionary, here we drop them into the proper cells of the 2D arrays
    for k in range(bergMaski):
        j = bergDict[k+1][0]
        i = bergDict[k+1][1]
        bergDepths[:,j,i] = icebergs_depths[k,:]
        bergWidths[:,j,i] = icebergs_widths[k,:]
        bergLength[:,j,i] = icebergs_length[k,:]
        bergsPerCell[j,i] = icebergs_per_cell[k]

print('\t\tMelange edge at i=%i' %melangeIndex)
for i in range(nx): #loops over whole domain, this calculates new hfill factors
    for j in range(ny):
        numberOfBergs = len(bergDepths[bergDepths[:,j,i] > 0,j,i])
        # print(i,j,numberOfBergs)
        if(numberOfBergs > 0):
            if(i < (melangeIndex - 1)): #Only scale for bergs within melange, not in shadedEdge
                scaleFactor =  thicknessLookup(x[0,i]) * iceDensity/oceanDensity * icebergCoverLambda / effectiveDepth[i]
                scaleFactor = (scaleFactor if abs(scaleFactor - 1) > scaleMin else 1)
                if(scaleFactor != 1 and j == int(ny/2)):
                    print('\t\t\tbergs are adjusted: %.3f i=%i, nB: %i' %(scaleFactor,i,np.sum(bergsPerCell[:,i])))
                bergDepths[:,j,i] = bergDepths[:,j,i] * scaleFactor
                bergDepths[bergDepths[:,j,i] > maxBergDepth,j,i] = maxBergDepth - random.random() # set max depth, bump em down a bit. Integer depth can cause problems
                bergDepths[(bergDepths[:,j,i] > 0) & (bergDepths[:,j,i] < minBergDepth),j,i] = minBergDepth + random.random() # set min depth, bump em up a bit can cause problems
                bergDepths[bergDepths[:,j,i] > hardMaxDepth,j,i] = hardMaxDepth # this ensures we don't have crashing issues
            lengths = bergLength[bergLength[:,j,i] > 0,j,i] #return only non-zeros 
            widths = bergWidths[bergWidths[:,j,i] > 0,j,i] #return only non-zeros
            depths = bergDepths[bergDepths[:,j,i] > 0,j,i] #return only non-zeros (negatives should have been blocked before here)
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
            if(i > melangeIndex): #free drifting bergs
                barrierMask[j,i] = 0
                driftMask[j,i] = 1
            else: # in mélange
                barrierMask[j,i] = 1
                driftMask[j,i] = 0
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


#remove blocking from plume locations
barrierMask[plumeMask[:,1] > 1,1] = 0

# Visual outputs, can toggle on/off or save

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
    plt.plot(x[0,:],-effectiveDepth/icebergCoverLambda,color='xkcd:red',linestyle='--',label='$D_{eff}$')
    plt.plot(x[0,:],-np.nanmean(depthHelper[:,1:-1,:],axis=(0,1)),color='xkcd:orange',linestyle='--',label='$<D>$')
plt.plot(X_,-H*iceDensity/oceanDensity,color='xkcd:green',linestyle='--',label='GLACIOME profile')
plt.plot([glaciome_L,glaciome_L],[20,-600],linestyle='--',color='cyan',label='Mélange Front')
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
# plt.savefig('figs/advectBergs%08d.png' %maxStep,format='png',dpi=150)

# I toggle this for troubleshooting purposes
if(not trialAndPlot):
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
