#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import numpy as np

# Input Variables here appropriately
nx = 50
ny = 12
nz = 20

dz = 10 # [m]
dx = 500
dy = 500

# How do you want to adjust the Bergs?
scale = .7
minThick = 12.00


# In[ ]:


#2D Fields
bergMask = np.fromfile('bergMask.bin', dtype='>f8')
bergMask = bergMask.reshape((ny,nx))
bergMaskNums = np.fromfile('bergMaskNums.bin', dtype='>f8')
bergMaskNums = bergMaskNums.reshape((ny,nx))
bergsPerCell = np.fromfile('numBergsPerCell.bin', dtype='>f8')
bergsPerCell = bergsPerCell.reshape((ny,nx))
driftMask = np.fromfile('driftMask.bin', dtype='>f8')
driftMask = driftMask.reshape((ny,nx))
barrierMask = np.fromfile('barrierMask.bin', dtype='>f8')
barrierMask = barrierMask.reshape((ny,nx))
#3D Fields
totalBergArea = np.fromfile('totalBergArea.bin', dtype='>f8')
totalBergArea = totalBergArea.reshape((nz,ny,nx))
openFrac = np.fromfile('openFrac.bin', dtype='>f8')
openFrac = openFrac.reshape((nz,ny,nx))

os.rename("totalBergArea.bin", "totalBergArea.bin" + "~") #save old file incase something goes wrong
os.rename("openFrac.bin", "openFrac.bin" + "~") #save old file incase something goes wrong


# In[ ]:


#Adjust Depth of the bergs
for i in range(nx):
    for j in range(ny):
        bergCount = int(bergsPerCell[j,i])
        if(bergMask[j,i] == 1 and bergCount > 0):  #only go in if bergs here
            #print('Berg number', bergMaskNums[j,i])
            depthFile = 'iceberg_depth_%05i.txt' % int(bergMaskNums[j,i])
            # widthFile = 'iceberg_width_%05i.txt' % int(bergMaskNums[j,i])
            # lengthFile = 'iceberg_length_%05i.txt' % int(bergMaskNums[j,i])
            os.rename(depthFile, depthFile + "~") #save old file incase something goes wrong
            with open(depthFile + "~",'r') as readFile:
                with open(depthFile, "w") as writeFile:
                    for line in readFile:
                        writeFile.write(str(max(float(line)*scale,minThick)) + "\n")
            readFile.close()
            writeFile.close()
            os.remove(depthFile + "~")
            


# In[ ]:


# Now we went and messed up the bergs, we need to adjust the global values
del depthFile

openFracTest = np.zeros([nz,ny,nx])
SATest = np.zeros([nz,ny,nx])
SATest[:,:,:] = np.nan
cellVolume = dz*dx*dy
for i in range(nx):
    for j in range(ny):
        bergCount = int(bergsPerCell[j,i])
        if(bergMask[j,i] == 1 and bergCount > 0):  #only go in if bergs here
            #print('Berg number', bergMaskNums[j,i])
            depthFile = 'iceberg_depth_%05i.txt' % int(bergMaskNums[j,i])
            widthFile = 'iceberg_width_%05i.txt' % int(bergMaskNums[j,i])
            lengthFile = 'iceberg_length_%05i.txt' % int(bergMaskNums[j,i])
            
            depths = np.zeros(bergCount)
            widths = np.zeros(bergCount)
            lengths = np.zeros(bergCount)
            #print('Berg count', bergCount)
            with open(depthFile,'r') as readFile:
                ii = 0
                for line in readFile:
                    if ii >= bergCount:
                        print('berg count mismatch in depth')
                        break
                    depths[ii] = float(line)
                    ii += 1
            readFile.close()
            with open(widthFile,'r') as readFile:
                ii = 0
                for line in readFile:
                    if ii >= bergCount:
                        print('berg count mismatch in width')
                        break
                    widths[ii] = float(line)
                    ii += 1
            readFile.close()
            with open(lengthFile,'r') as readFile:
                ii = 0
                for line in readFile:
                    if ii >= bergCount:
                        print('berg count mismatch in length')
                        break
                    lengths[ii] = float(line)
                    ii += 1
            readFile.close()
            #print(depths,lengths,widths)
            for k in range(nz):
                d_bot = k*dz + dz #bottom of depth bin
                d_top = k*dz
                volume1 = dz * lengths[depths > d_bot] * widths[depths > d_bot]
                SA1 = dz*2*(lengths[depths > d_bot] + widths[depths > d_bot])
                partialFill = (depths < d_bot) & (depths > d_top)
                #partial fill
                volume2 = (depths[partialFill] - d_top) * lengths[partialFill] * widths[partialFill]
                #partial sides
                SA2 = (depths[partialFill] - d_top)*2*(lengths[partialFill] + widths[partialFill]) 
                #bottom
                SA3 = lengths[partialFill] * widths[partialFill]
                #print(np.sum(volume1), np.sum(volume2))
                openFracTest[k,j,i] = 1-((np.sum(volume1) + np.sum(volume2))/cellVolume)
                SATest[k,j,i] = np.sum(SA1) + np.sum(SA2) + np.sum(SA3)
        elif(bergMask[j,i] == 1 and bergCount == 0):
            openFracTest[:,j,i] = 1
            SATest[:,j,i] = 0

print('\\( ﾟヮﾟ)/ Done! \nAll your icebergs are now %i%% the size!' % (scale*100))
print('openFrac.bin and totalBergArea.bin have been updated to reflect the change.')


# In[ ]:


writeFiles = True
def write_bin(fname, data):
    print(fname, np.shape(data))
    
    if(writeFiles):
        data.astype(">f8").tofile(fname)
    else:
        setupNotes.write('Not saving')


# swaped = np.swapaxes(np.swapaxes(openFracTest,0,0),0,0)
write_bin('openFrac.bin',openFracTest)
write_bin('totalBergArea.bin',SATest)


# In[ ]:


openFrac = np.fromfile('openFrac.bin~', dtype='>f8')
openFracShaped = openFrac.reshape((nz,ny,nx))
openFracN = np.fromfile('openFrac.bin', dtype='>f8')
openFracNShaped = openFracN.reshape((nz,ny,nx))
SA = np.fromfile('totalBergArea.bin~', dtype='>f8')
SAShaped = openFrac.reshape((nz,ny,nx))
SAN = np.fromfile('totalBergArea.bin', dtype='>f8')
SANShaped = openFracN.reshape((nz,ny,nx))


# In[ ]:


print('Max volume error:',np.max((openFracNShaped[:,:,0:8]-openFracShaped[:,:,0:8])/(openFracShaped[:,:,0:8]+ 1e-8)))
print('Max SA error:',np.nanmax((SATest[:,:,0:8]-totalBergArea[:,:,0:8])/(totalBergArea[:,:,0:8]+ 1e-8)))


# In[ ]:





# In[ ]:




