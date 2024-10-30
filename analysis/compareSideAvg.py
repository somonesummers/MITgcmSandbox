from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import fileinput
import cmocean

# Pick cross section to view
crossSection = 1000
try:
    with open('input/plotPoint.txt', 'r') as file:
        lines = file.readlines()
        crossSection = float(lines[2]) #reads the 3rd line in the doc
        print('cross section read from file', crossSection)
except FileNotFoundError:
    print('plot point file does not exist, using default')

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

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


if((maxStep-startStep)/sizeStep > 60):   #if more than 50 frames, downscale to be less than 50
    dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/60)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

#Decide if iceBerg data files exist
if(os.path.isfile('input/bergMask.bin')):
    isBerg = True
    print('Found icebergs for this run')
else:
    isBerg = False

os.system('rm -f figs/compareSide_*.png')
os.system('rm -f figs/compareSide_*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")

if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))


if(isBerg):
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    bergMaskNums = np.fromfile('input/bergMaskNums.bin', dtype='>f8')
    bergMaskNums = bergMaskNums.reshape(np.shape(x))
    bergsPerCell = np.fromfile('input/numBergsPerCell.bin', dtype='>f8')
    bergsPerCell = bergsPerCell.reshape(np.shape(x))
    bergContaingingCells = int(np.sum(bergMask))
    maxDepth = np.zeros(np.shape(x))

    #deepest contour
    for i in range(np.shape(x)[1]):
        for j in range(np.shape(x)[0]):
            bergCount = int(bergsPerCell[j,i])
            if(bergMask[j,i] == 1 and bergCount > 0):  #only go in if bergs here
                depthFile = 'input/iceberg_depth_%05i.txt' % int(bergMaskNums[j,i])
                depths = np.zeros(bergCount)
                with open(depthFile,'r') as readFile:
                    ii = 0
                    for line in readFile:
                        if ii >= bergCount:
                            print('berg count mismatch in depth')
                            break
                        depths[ii] = float(line)
                        ii += 1
                readFile.close()
                maxDepth[j,i] = np.max(depths)
    # contourf plot
    openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1
                    

ySlice = np.argmin(np.abs(y[:,0] - crossSection))
print('cross section is y =', y[ySlice,0], 'index', ySlice)

dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag']
name = ["Temp", "Sal", "U", "W", "V"]
cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]"]

for k in range(len(name)):
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        if(isBerg and os.path.isfile('results/BRGFlx.%010i.001.001.data' % i)):
            localBergs = True
        else:
            localBergs = False
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        data_old = mds.rdmds("results/%s"%(dynName[k]), 20*86400/dt)
        data = data-data_old
        if k == 0:  #T
            lvl = np.linspace(-1.5, 1.5, 127)
            cm = "cmo.tarn_r"
        elif k == 1: #S
            lvl = np.linspace(-0.5, 0.5, 127)
            cm = "cmo.diff"
        elif k == 2 or k == 4: #u,v
            lvl = np.linspace(-.1, .1, 127)
            cm = "cmo.balance"
        elif k == 3: #W
            lvl = np.linspace(-0.01, 0.01, 127)
            cm = "cmo.curl"
        elif k == 5:
            lvl = np.linspace(0, .5, 128)
            cm = "cmo.rain"
        if(k == 5):
            kk = 2
        else:
            kk = k
        cp = plt.contourf(
            np.squeeze(x[ySlice,:]),
            np.squeeze(z),
            np.squeeze(data[kk, :, ySlice, :]),
            lvl,
            extend="both",
            cmap=cm,
        )
        plt.plot(x[ySlice,:],topo[ySlice,:],color='black')
        if(localBergs):
            plt.plot(x[ySlice,:],-np.max(maxDepth,axis=0),color='gray',linestyle='dotted')
            cp2 = plt.contourf(
                x[ySlice,:],
                np.squeeze(z),
                np.squeeze(openFrac[:, ySlice, :]),
                [.1,.5,.9],
                extend="min",
                alpha=.1,
                cmap='cmo.gray')
            # cbar2 = plt.colorbar(cp2)
            # cbar2.set_label('Ocean Fraction')
        cbar = plt.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(data[kk, :, ySlice, :]),np.nanmax(data[kk, :, ySlice, :]),np.max(np.isnan(data[kk, :, ySlice, :]))))
        plt.ylabel('Depth [m]')
        plt.title("%s y = %i at %.02f days" % (name[k], y[ySlice,0], i/86400*dt))
        j = i/sizeStep + startStep
        
        str = "figs/compareSide_%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png')
        plt.close()
        #plt.show()

    os.system('magick -delay %f figs/compareSide_%s*.png -colors 256 -depth 256 figs/compareSide_%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

#clean up intermediate pngs
os.system('rm -f figs/compareSide_*.png')
