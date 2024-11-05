from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import cmocean
import argparse

parser = argparse.ArgumentParser(description='Compare 2 MITcgm folders.')
parser.add_argument('directories', metavar='Dir', type=str, nargs=2,
                    help='2 folders to compare')
args = parser.parse_args()

#folder1 controls plotting steps, any bergs, bathymetry, etc
folder1 = args.directories[0]
folder2 = args.directories[1]

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = 50
plotDPI = 100

if(os.path.isfile('input/plotHelper.py')):
    sys.path.append('input')
    from plotHelper import *
    print('Found experiment plotting settings')
elif(os.path.isfile('plotHelper.py')):
    print('no custom plotting settings, using local default')
    from plotHelper import *
else:  
    print('no defaults found')

maxStep = 0
sizeStep = 1e10
startStep = 1e10

for file in os.listdir('%s/results' %folder1):
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


if(maxStep/sizeStep > 50):   #if more than 50 frames, downscale to be less than 50
    dwnScale = round((maxStep/sizeStep)/50)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)
#Decide if iceBerg data files exist
if(os.path.isfile('%s/input/data.iceberg' % folder1)):
    isBerg = True
else:
    isBerg = False

os.system('rm -f %s/figs/compareSide_*.png' % folder1)
os.system('rm -f %s/figs/autoCompareSide_*.gif' % folder1)

x = mds.rdmds("%s/results/XC" % folder1)
y = mds.rdmds("%s/results/YC" % folder1)
z = mds.rdmds("%s/results/RC" % folder1)

if(os.path.isfile('%s/input/bathymetry.bin' % folder1)):
    topo = np.fromfile('%s/input/bathymetry.bin' % folder1, dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))
# ice = np.fromfile('input/icetopo.exp1', dtype='>f8')

# ice = ice.reshape(np.shape(x))

if(isBerg):
    bergMask = np.fromfile('%s/input/bergMask.bin' % folder1, dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    bergMaskNums = np.fromfile('%s/input/bergMaskNums.bin' % folder1, dtype='>f8')
    bergMaskNums = bergMaskNums.reshape(np.shape(x))
    bergsPerCell = np.fromfile('%s/input/numBergsPerCell.bin' % folder1, dtype='>f8')
    bergsPerCell = bergsPerCell.reshape(np.shape(x))
    bergContaingingCells = int(np.sum(bergMask))
    maxDepth = np.zeros(np.shape(x))

    #deepest contour
    for i in range(np.shape(x)[1]):
        for j in range(np.shape(x)[0]):
            bergCount = int(bergsPerCell[j,i])
            if(bergMask[j,i] == 1 and bergCount > 0):  #only go in if bergs here
                depthFile = '%s/input/iceberg_depth_%05i.txt' % (folder1, int(bergMaskNums[j,i]))
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
    openFrac = np.fromfile('%s/input/openFrac.bin' % folder1, dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile('%s/input/bergMask.bin' % folder1, dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1
                    

ySlice = np.argmin(np.abs(y[:,0] - crossSection))
print('cross section is y =', y[ySlice,0], 'index', ySlice)

# if(isBerg):
#     dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'BRGFlx']
#     name = ["Temp", "Sal", "U", "W", "V","BRGmltRt"]
#     cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[m/d]"]
# else:
dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag']
name = ["Temp", "Sal", "U", "W", "V"]
cbarLabel = ["[∆ C]", "[∆ ppt]", "[∆ m/s]", "[∆ m/s]", "[∆ m/s]"]

for k in range(len(name)):
    #print('k,',k)
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data1 = mds.rdmds("%s/results/%s" % (folder1, dynName[k]), i)
        data2 = mds.rdmds("%s/results/%s" % (folder2, dynName[k]), i)
        data = data1 - data2
        if k == 0:  #T
            lvl = np.linspace(-1.5, 1.5, 128)
            cm = "cmo.tarn_r"
        elif k == 1: #S
            lvl = np.linspace(-0.5, 0.5, 128)
            cm = "cmo.diff"
        elif k == 2 or k == 4: #u,v
            lvl = np.linspace(-.1, .1, 128)
            cm = "cmo.balance"
        elif k == 3: #W
            lvl = np.linspace(-0.01, 0.01, 128)
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
        
        str = "%s/figs/compareSide_%s%05i.png" % (folder1, name[k],j)
        
        plt.savefig(str, format='png')
        plt.close()
        #plt.show()

    os.system('magick -delay %f %s/figs/compareSide_%s*.png -colors 256 -depth 256 %s/figs/autoComapreSide_%s.gif' %(500/(maxStep/sizeStep), folder1, name[k], folder1, name[k]))

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
