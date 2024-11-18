#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os
from matplotlib import pyplot as plt
import time
import cmocean
writeFiles = True

# Cell resolution
deltaX = 500
deltaY = 500
deltaZ = 20

# Dimensions of grid
nx = 80 # 15 km long
ny = 12 # 10 km wide (plus 500 m walls)
nz = 50 # 1000 m deep

delz = np.ones(nz)*deltaZ

#Make Masks
bergMask = np.zeros([ny,nx])
driftMask = np.zeros([ny,nx])
barrierMask = np.zeros([ny,nx])
bergConc = np.zeros([ny,nx])

bergMaskNums = np.zeros([ny,nx])
numBergsPerCell = np.zeros([ny,nx],dtype=np.int16)

# Berg parameters
bergType = 1 # 1 = block 2 = cone (not implemented)
alpha = 1.8 # slope of inverse power law size frequency distribution
scaling = 2 # 1 = Sulak 2017 2 = Barker 2004
maxBergDepth = 300 # (m) - set to zero if 'prescribing' max iceberg width
minBergDepth= 20 # (m)
maxBergWidth = 0 # (m) - set to zero if 'prescribing' max iceberg depth
minBergWidth = 20 # (m)

# Iceberg mask
bergMask[1:-1,1:20] = 1 # icebergs in inner 5 km, all oriented east-west

# Drift mask
driftMask[1:-1,1:20] = 1 # calculate effect of iceberg drift on melt rates 

# Barrier mask
barrierMask[1:-1,1:20] = 1 # make icebergs a physical barrier to water flow

# Iceberg concentration (# of each surface cell that is filled in plan view)
bergConc[1:-1,1:20] = 95 # iceberg concentration 

desiredBergArea = np.sum(bergConc/100.0*deltaX*deltaY)
bergMaskArea = np.sum(bergMask*deltaX*deltaY)
print('Area where bergs live:',bergMaskArea,'m^2')
print('Desired berg area:',desiredBergArea,'m^2')
print('Ratio:',desiredBergArea/bergMaskArea*100,'%')


# In[3]:


if(scaling == 1): # then use Sulak17 volume-area scaling volume = 6.0*area^1.30
    # assumes volume = L*W*D and W = L/1.62 (Dowdeswell et al 1992)
    if(maxBergWidth==0):
        maxBergWidth = 0.0642449*maxBergDepth**(5/3)
        # minBergWidth = 0.0642449*minBergDepth**(5/3)
    elif(maxBergDepth==0):
        maxBergDepth = 5.19155*maxBergWidth**(5/3)
        minBergDepth = 5.19155*minBergWidth**(5/3)
elif(scaling == 2): # Then use Barker04 width-depth relationship
    # Depth = 2.91*Width^0.71
    if(maxBergWidth==0):
        maxBergWidth = (100*10**(58/71)*maxBergDepth**(100/71)) / (291*291**(29/71))
        #minBergWidth = (100*10**(58/71)*minBergDepth**(100/71)) / (291*291**(29/71))        
    elif(maxBergDepth==0):
        maxBergDepth = 2.91*maxBergWidth^0.71
        minBergDepth = 2.91*minBergWidth^0.71


# In[4]:


numberOfBergs = 1000 #for now, immediately doubled by scheme below, so guess low, high guesses (300%+) can cause to fail
bergTopArea = 0
areaResidual = 1
# Generate the Inverse Power Law cumulative distribution function
# over the range minBergWidth-maxBergWidth with a slope of alpha.

while(np.abs(areaResidual) > .01 ):

    numberOfBergs = round(numberOfBergs * (1 + areaResidual))
    print('numberOfBergs',numberOfBergs)
    x_width = np.arange(minBergWidth, maxBergWidth, (maxBergWidth-minBergWidth)/(numberOfBergs*1e2))
    x_depth = np.arange(minBergDepth, maxBergDepth, (maxBergDepth-minBergDepth)/(numberOfBergs*1e2))
    inversePowerLawPDF_width = ((alpha-1) / minBergWidth) * (x_width/minBergWidth) ** (-alpha)
    inversePowerLawPDF_depth = ((alpha-1) / minBergDepth) * (x_depth/minBergDepth) ** (-alpha)
        # Get the CDF numerically
    inversePowerLawCDF_width = np.cumsum(inversePowerLawPDF_width)
    inversePowerLawCDF_depth = np.cumsum(inversePowerLawPDF_depth)
        # Normalize
    inversePowerLawCDF_width = inversePowerLawCDF_width / inversePowerLawCDF_width[-1]
    inversePowerLawCDF_depth = inversePowerLawCDF_depth / inversePowerLawCDF_depth[-1]
        
        # Generate number_of_bergs uniformly distributed random numbers.
    uniformlyDistributedRandomNumbers = np.random.uniform(0,1,numberOfBergs)
    
    # plt.subplot(1,2,1)
    # plt.plot(inversePowerLawCDF_width)
    
    # plt.subplot(1,2,2)
    # plt.plot(inversePowerLawCDF_depth)
    # Invert the CDF of the Inverse Power Law function to get a function that can
        # generate random numbers drawn from a Inverse Power Law distribution,
        # given numbers drawn from a uniform distribution.
    
    inversePowerLawDistNumbers_width = np.zeros(uniformlyDistributedRandomNumbers.size);
    inversePowerLawDistNumbers_depth = np.zeros(uniformlyDistributedRandomNumbers.size);
    nearestIndex_width = [0] * uniformlyDistributedRandomNumbers.size
    nearestIndex_depth = [0] * uniformlyDistributedRandomNumbers.size
    
    for i in range(uniformlyDistributedRandomNumbers.size):  #this is pretty slow, not sure if vectorized exists
        nearestIndex_width[i] = np.abs(uniformlyDistributedRandomNumbers[i]-inversePowerLawCDF_width).argmin();
        nearestIndex_depth[i] = np.abs(uniformlyDistributedRandomNumbers[i]-inversePowerLawCDF_depth).argmin();
    
    inversePowerLawDistNumbers_width = x_width[nearestIndex_width];
    inversePowerLawDistNumbers_length = inversePowerLawDistNumbers_width/1.62 # Widths are bigger 
    tooWide = np.count_nonzero(inversePowerLawDistNumbers_width > deltaX)
    tooLong = np.count_nonzero(inversePowerLawDistNumbers_width > deltaX)
    inversePowerLawDistNumbers_width[inversePowerLawDistNumbers_width > deltaX] = deltaX # Max width is grid cell (assumed square)
    inversePowerLawDistNumbers_length[inversePowerLawDistNumbers_length > deltaX] = deltaX # Max length is grid cell (assumed square)
    if(tooLong + tooWide > 0):
        print('Bergs clipped: %i for width, %i for length' % (tooWide, tooLong))
    
    inversePowerLawDistNumbers_depth = x_depth[nearestIndex_depth]; #depths don't get clipped
    
    bergTopArea = sum(inversePowerLawDistNumbers_width*inversePowerLawDistNumbers_length)
    areaResidual = (desiredBergArea - bergTopArea)/desiredBergArea
    print(bergTopArea/bergMaskArea*100,"% Bergs")
    print('areaResidual (%):', areaResidual)

print('====== Success! Found our bergs =====')
print('Width min/mean/max:',np.min(inversePowerLawDistNumbers_width),np.mean(inversePowerLawDistNumbers_width),np.max(inversePowerLawDistNumbers_width),'[meters]')
print('Depth min/mean/max:',np.min(inversePowerLawDistNumbers_depth),np.max(inversePowerLawDistNumbers_depth),np.max(inversePowerLawDistNumbers_depth),'[meters]')
print('Total Berg Area', bergTopArea)
print('Total Berg fract:', bergTopArea/bergMaskArea,'%')


# In[5]:


bergMaski = 0  #Bad name, but this is the count of cells that will recieve bergs
bergDict = {}

for j in range(ny):
    for i in range(nx):
        if(bergMask[j,i] == 1):
            # print('i,j, bergmask',i,j,bergMask[j,i])
            bergMaski = 1 + bergMaski #Needs to start at 1, as non-bergs will be 0
            bergMaskNums[j,i] = bergMaski #Assign Mask Nums, not random as we'll randomly place bergs in cells
            bergDict[bergMaski] = [j,i] #This lets us do 1-D loops for the whole grid

print(bergMaski, 'cells with bergs')
# print(bergDict)

# Sort my bergs
sorted_indices = np.argsort(-inversePowerLawDistNumbers_depth) # Sort backwards to get descending from big to small bergs 
sorted_depth = inversePowerLawDistNumbers_depth[sorted_indices]
sorted_width = inversePowerLawDistNumbers_width[sorted_indices]
sorted_length = inversePowerLawDistNumbers_length[sorted_indices]

# print(np.sum(inversePowerLawDistNumbers_width * inversePowerLawDistNumbers_length))
# print(np.sum(sorted_width * sorted_length))

assignedCell = np.random.randint(0,bergMaski,[numberOfBergs]) # In this script, every berg has a home


# In[35]:


# Need to reshuffle as it re-writes this when sorting
# Array for bergs
bergsPerCellLimit = 5000
icebergs_depths = np.zeros([bergMaski,bergsPerCellLimit])
icebergs_widths = np.zeros([bergMaski,bergsPerCellLimit])
icebergs_length = np.zeros([bergMaski,bergsPerCellLimit])

assignedCell = np.random.randint(0,bergMaski,[numberOfBergs]) #every Berg has a spot
print('max cell assigned bergs', np.max(assignedCell))

icebergs_per_cell = np.zeros([bergMaski],dtype=np.int16)
icebergs_area_per_cell = np.zeros([bergMaski])
# numberOfBergs = sorted_indices.size  this is now taken care of appropriately above

for i in range(numberOfBergs): 
    j = assignedCell[i]
    # print('looking at mask number',j,'at berg',i)
    # print('Berg number', icebergs_per_cell[j],'in this cell')
    bergArea = sorted_width[i] * sorted_length[i]
    loopLimiter = 0
    while(bergArea > (deltaX * deltaY  * (bergConc[bergDict[j+1][0],bergDict[j+1][1]]/100 + .02) - icebergs_area_per_cell[j])): #if too full, pick random new cell
        j_old = j
        j = np.random.randint(0,bergMaski)
        assignedCell[i] = j
        loopLimiter += 1
        # print(j_old, "is full, trying",j, "attempt", loopLimiter)
        if(loopLimiter > bergMaski):
            print('no more empty cells... giving up')
            break
    icebergs_depths[j,icebergs_per_cell[j]] = sorted_depth[i]
    icebergs_widths[j,icebergs_per_cell[j]] = sorted_width[i]
    icebergs_length[j,icebergs_per_cell[j]] = sorted_length[i]
    icebergs_per_cell[j] += 1
    # icebergs_area_per_cell[j] = np.sum(icebergs_widths[j,:]*icebergs_length[j,:])
    icebergs_area_per_cell[j] += bergArea
print('Bergs per cell and filled faction at surface for spot check')     
print(icebergs_per_cell)
print(np.round(icebergs_area_per_cell/(500*500),2))
# print('sum of areas:',np.sum(icebergs_area_per_cell))
# print('Other  areas:',np.sum(icebergs_widths * icebergs_length))
# All bergs now sorted 


# In[7]:


openFrac = np.zeros([nz,ny,nx])
SA = np.zeros([nz,ny,nx])
SA[:,:,:] = np.nan
cellVolume = deltaX*deltaY*deltaZ

#This loop knows where all bergs are already, different from searching for all bergs across entire grid
for i in range(bergMaski):
    bergCount = icebergs_per_cell[i]
    numBergsPerCell[bergDict[i+1][0],bergDict[i+1][1]] = bergCount
    if(bergCount > 1):
        lengths = icebergs_length[i,icebergs_length[i,:] > 0] #return only non-zeros
        widths = icebergs_widths[i,icebergs_widths[i,:] > 0] #return only non-zeros
        depths = icebergs_depths[i,icebergs_depths[i,:] > 0] #return only non-zeros
        for k in range(nz):
            d_bot = k*deltaZ + deltaZ #bottom of depth bin
            d_top = k*deltaZ
            volume1 = deltaZ * lengths[depths > d_bot] * widths[depths > d_bot]
            SA1 = deltaZ*2*(lengths[depths > d_bot] + widths[depths > d_bot])
            partialFill = (depths < d_bot) & (depths > d_top)
            #partial fill
            volume2 = (depths[partialFill] - d_top) * lengths[partialFill] * widths[partialFill]
            #partial sides
            SA2 = (depths[partialFill] - d_top)*2*(lengths[partialFill] + widths[partialFill]) 
            #bottom
            SA3 = lengths[partialFill] * widths[partialFill]
            #print(np.sum(volume1), np.sum(volume2))
            openFrac[k,bergDict[i+1][0],bergDict[i+1][1]] = 1-((np.sum(volume1) + np.sum(volume2))/cellVolume)
            SA[k,bergDict[i+1][0],bergDict[i+1][1]] = np.sum(SA1) + np.sum(SA2) + np.sum(SA3)
    elif(bergCount == 0):
        openFrac[:,bergDict[i+1][0],bergDict[i+1][1]] = 1
        SA[:,bergDict[i+1][0],bergDict[i+1][1]] = 0
# print(np.round(openFrac[0,:,0:8],2))


# In[34]:


# print('Number of bergs to cells for those curious')
# print(numBergsPerCell)


# In[9]:


fig = plt.figure()
plt.subplot(2,2,1)
for i in range(bergMaski):
    plt.plot(openFrac[:,bergDict[i+1][0],bergDict[i+1][1]],-np.cumsum(delz),alpha=.5,color='xkcd:gray',linewidth=.5)
plt.plot(np.mean(openFrac[:,bergMask==1],1),-np.cumsum(delz),alpha=1,color='xkcd:black',linewidth=1,linestyle='--',label='Average Bergs')
plt.plot([0,1],[-maxBergDepth,-maxBergDepth],color = 'xkcd:red',linestyle=':', label='Target Max Depth')
plt.plot([1-np.max(bergConc)/100,1-np.max(bergConc)/100],[-nz*deltaZ,0],color = 'xkcd:gray',linestyle=':',label='Target Max Berg Conc')
# plt.title('Open fraction of cells as function of depth across all berg containing Cells')
plt.xlabel('Open Fraction of Cells')
plt.ylabel('Depth [m]')
# plt.legend()

plt.subplot(2,2,2)
plt.hist(inversePowerLawDistNumbers_depth,bins = 50)
# plt.ylabel('Count')
plt.xlabel('Depth [m]')

plt.subplot(2,2,3)
plt.hist(inversePowerLawDistNumbers_width,bins = 50)
plt.ylabel('Count')
plt.xlabel('Width [m]')

plt.subplot(2,2,4)
plt.hist(inversePowerLawDistNumbers_length,bins = 50)
# plt.ylabel('Count')
plt.xlabel('Length [m]')
fig.tight_layout()
plt.savefig('bergStatistics.png', format='png', dpi=200)
plt.show()


# In[30]:


pltHelper = 1-openFrac[0,:,:]
pltHelper[bergMask == 0] = np.nan
fig = plt.figure()
plt.subplot(211)

pc = plt.pcolor(pltHelper,cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
plt.suptitle('Iceberg Cover and Residual')
plt.ylabel('Cell across fjord')
cbar.set_label('iceberg cover')

plt.subplot(212)
pc = plt.pcolor(pltHelper - bergConc/100,cmap='cmo.ice')
cbar = plt.colorbar(pc)

pc_min = str(np.round(np.nanmin(pltHelper),3) * 100)
pc_max = str(np.round(np.nanmax(pltHelper),3) * 100)
plt.xlabel("Cell along fjord, Coverage is " + pc_min + "%-" + pc_max + "%")
plt.ylabel('Cell across fjord')

cbar.set_label('cover resid')

fig.tight_layout()
plt.savefig('bergMap.png', format='png', dpi=200)
plt.show()


# In[11]:


def write_bin(fname, data):
    print(fname, np.shape(data))
    
    if(writeFiles):
        data.astype(">f8").tofile(fname)
    else:
        print('Not saving')

# write iceberg txt files
if(writeFiles):
    for i in range(bergMaski):
        with open('iceberg_depth_%05i.txt' % (i+1) , 'w') as file_handler:
            for item in icebergs_depths[i,icebergs_depths[i,:] > 0]:
                file_handler.write("{}\n".format(item))
        with open('iceberg_width_%05i.txt' % (i+1) , 'w') as file_handler:
            for item in icebergs_widths[i,icebergs_widths[i,:] > 0]:
                file_handler.write("{}\n".format(item))
        with open('iceberg_length_%05i.txt' % (i+1) , 'w') as file_handler:
            for item in icebergs_length[i,icebergs_length[i,:] > 0]:
                file_handler.write("{}\n".format(item))

# write global files
write_bin('bergMask.bin',bergMask)
write_bin('bergMaskNums.bin',bergMaskNums)
write_bin('numBergsPerCell.bin',numBergsPerCell)
write_bin('openFrac.bin',openFrac)
write_bin('totalBergArea.bin',SA)
write_bin('driftMask.bin',driftMask)
write_bin('barrierMask.bin',barrierMask)

print('Done making bergs.')