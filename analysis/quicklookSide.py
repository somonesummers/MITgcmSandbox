from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import cmocean

# Pick cross section to view
crossSection = 1000
with open('input/plotPoint.txt', 'r') as file:
    lines = file.readlines()
    crossSection = float(lines[2]) #reads the 3rd line in the doc
    print('cross section read from file', crossSection)


maxStep = 0
sizeStep = 1e10
startStep = 1e10

for file in os.listdir('results'):
    # print(file)
    if "dynDiag" in file:
        words = file.split(".")
        # print(words[1])  
        if int(words[1]) > maxStep:
            maxStep = int(words[1])
        if int(words[1]) < startStep and int(words[1]) > 0:
            sizeStep = int(words[1])
            startStep = int(words[1])

if(maxStep/sizeStep > 50):  #if more than 50 frames, downscale to be less than 50
    dwnScale = round((maxStep/sizeStep)/50)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print(startStep,sizeStep,maxStep)

#Decide if iceBerg data files exist
if(os.path.isfile('input/data.iceberg')):
    isBerg = True
else:
    isBerg = False

os.system('rm -f figs/side_*.png')
os.system('rm -f figs/autoside_*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")

if(os.path.isfile('input/topog.slope')):
    topo = np.fromfile('input/topog.slope', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))
# ice = np.fromfile('input/icetopo.exp1', dtype='>f8')

# ice = ice.reshape(np.shape(x))

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

if(isBerg):
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'BRGFlx']
    name = ["Temp", "Sal", "U", "W", "V","BRGmltRt"]
    cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[m/d]"]
else:
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag']
    name = ["Temp", "Sal", "U", "W", "V"]
    cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]"]

for k in range(len(name)):
    #print('k,',k)
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        if k == 0:
            lvl = np.linspace(-1.5, 1.5, 128)
            cm = "cmo.thermal"
        elif k == 1:
            lvl = np.linspace(32, 34, 128)
            cm = "cmo.haline"
        elif k == 2 or k == 4:
            lvl = np.linspace(-.5, .5, 128)
            cm = "cmo.balance"
        elif k == 3:
            lvl = np.linspace(-0.08, 0.08, 128)
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
        # plt.plot(x[ySlice,:],ice[ySlice,:],color='gray')
        if(isBerg):
            plt.plot(x[ySlice,:],-np.max(maxDepth,axis=0),color='gray',linestyle='dotted')
            cp2 = plt.contourf(
                x[ySlice,:],
                np.squeeze(z),
                np.squeeze(openFrac[:, ySlice, :]),
                [.1,.5,.9],
                extend="min",
                alpha=.2,
                cmap='cmo.gray')
            cbar2 = plt.colorbar(cp2)
            cbar2.set_label('Ocean Fraction')
        cbar = plt.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(data[kk, :, ySlice, :]),np.nanmax(data[kk, :, ySlice, :]),np.max(np.isnan(data[kk, :, ySlice, :]))))
        plt.ylabel('Depth [m]')
        plt.title("%s y = %i at %i" % (name[k], y[ySlice,0], i))
        j = i/startStep
        
        str = "figs/side%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png')
        plt.close()
        #plt.show()

    os.system('magick -delay 5 figs/side%s*.png -colors 256 -depth 256 figs/autoside_%s.gif' %(name[k], name[k]))

# BCT = np.fromfile("T.bound", dtype=">f8")
# plt.plot(BCT)
# plt.show()
#
# BCS = np.fromfile("S.bound", dtype=">f8")
# plt.plot(BCS)
# plt.show()
#
# BCU = np.fromfile("U.bound", dtype=">f8")
# plt.plot(BCU)
# plt.show()
#
#
