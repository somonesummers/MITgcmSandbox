from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import argparse

parser = argparse.ArgumentParser(description='Compare 2 MITcgm folders, point to results folders.')
parser.add_argument('directories', metavar='Dir', type=str, nargs=2,
                    help='2 results folders to compare')
args = parser.parse_args()

#folder1 controls plotting steps, any bergs, bathymetry, etc
folder1 = args.directories[0]
folder2 = args.directories[1]

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 100
cleanPNGs = True
usePcolor = True
showQuiver = True
showZeros = True

if(os.path.isfile('input/plotHelperLocal.py')):
    sys.path.append('input')
    from plotHelperLocal import *
    print('Found experiment plotting settings')
elif(os.path.isfile('../plotHelper.py')):
    sys.path.append('../')
    print('no custom plotting settings, using local default')
    from plotHelper import *
else:  
    print('no defaults found')
print('Plot DPI:',plotDPI,'; clean PNGs:',cleanPNGs, '; usePcolor:', usePcolor)

dt1 = 0.0   
for line in fileinput.input('input/data'):
    if "deltaT=" in line:
        dt1 = float(line[8:-2])

dt2 = 0.0   
for line in fileinput.input(folder2 + '/../input/data'):
    if "deltaT=" in line:
        dt2 = float(line[8:-2])
print('dt 1 is loaded as', dt1,'dt 2 is loaded as', dt2)

if(dt1 != dt2):
    print("\tUneven timesteps detected, will attempt to find same timestep diagnostic files")

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
if(os.path.isfile(folder1 + '/../input/bergMask.bin')):
    isBerg1 = True
    print('Found icebergs for folder1')
else:
    isBerg1 = False
if(os.path.isfile(folder2 + '/../input/bergMask.bin')):
    isBerg2 = True
    print('Found icebergs for folder2')
else:
    isBerg2 = False

#Clean up old gifs and pngs
os.system('rm -f figs/compareMelt*.png')
os.system('rm -f figs/compareMelt*.gif')

y = mds.rdmds("results/YC")
x = mds.rdmds("results/XC")
z = mds.rdmds("results/RC")



if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))

if(isBerg1):
    # contourf plot
    openFrac = np.fromfile(folder1 + '/../input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile(folder1 + '/../input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1
elif(isBerg2):
    # contourf plot
    openFrac = np.fromfile(folder2 + '/../input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile(folder2 + '/../input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1


zSlice = np.argmin(np.abs(z[:,0,0]- zDepth))
print('depth is z =', z[zSlice,0,0], 'index', zSlice, 'for quiver plotting')


name = ["fwFlx", "htFlx", "Tau_x", "Tau_y"]
cbarLabel = ["[∆ m^3/s]", "[∆ W/m^2]", "[∆ Pa/m^2]", "[∆ Pa/m^2]"]

for k in range(len(name)):
    print("\t" + name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        if(showQuiver):
            dataQuiv = mds.rdmds(folder1 + "/dynDiag", i)
        if(isBerg1):
            if(k == 0):
                data1 = np.sum(mds.rdmds(folder1 + "/BRGFlx", i)[k,:,:,:],0)
            elif(k == 1):
                data1 = np.sum(mds.rdmds(folder1 + "/BRGFlx", i)[0,:,:,:],0) * 998 * 333.55E+3 / ((x[1,1] - x[1,0]) * (y[1,1] - y[0,1]))
            else:
                data1 = np.sum(mds.rdmds(folder1 + "/BRGFlx", i)[k+1,:,:,:],0)
        else:
             # raw data is kg/m^2/s, convert to m^3/s
            if(k == 0):
                data1 = mds.rdmds(folder1 + "/shelfDiag", i)[0,:,:] * -1 / 1000 * (x[1,1] - x[1,0]) * (y[1,1] - y[0,1])
            elif(k == 1):
                # raw data is W/m^2
                data1 = mds.rdmds(folder1 + "/shelfDiag", i)[k,:,:]
            else:
                data1 = mds.rdmds(folder1 + "/shelfDiag", i)[k,:,:]
        if(isBerg2):
            if(k == 0 ):
                data2 = np.sum(mds.rdmds(folder2 + "/BRGFlx", i * dt1 / dt2)[k,:,:,:],0)
            elif(k == 1):
                data2 = np.sum(mds.rdmds(folder2 + "/BRGFlx", i)[0,:,:,:],0) * 998 * 333.55E+3 / ((x[1,1] - x[1,0]) * (y[1,1] - y[0,1]))
            else:
                data2 = np.sum(mds.rdmds(folder2 + "/BRGFlx", i * dt1 / dt2)[k+1,:,:,:],0)
        else:
             
            if(k == 0):
                # raw data is kg/m^2/s, convert to m^3/s
                data2 = mds.rdmds(folder2 + "/shelfDiag", i * dt1 / dt2)[0,:,:] * -1 / 1000 * (x[1,1] - x[1,0]) * (y[1,1] - y[0,1])
            elif(k == 1):
                # raw data is W/m^2
                data2 = mds.rdmds(folder2 + "/shelfDiag", i * dt1 / dt2)[k,:,:]
            else:
                data2 =  mds.rdmds(folder2 + "/shelfDiag", i * dt1 / dt2)[k,:,:]    
        dataPlot = data1 - data2    
        if(k == 0):
            lvl = np.linspace(-2, 2, 127)
            cm = "cmo.tarn_r"
        elif(k == 1):
            lvl = np.linspace(-2500, 2500, 127)
            cm = "cmo.balance"
        elif(k == 2):
            lvl = np.linspace(-10000, 10000, 127)
            cm = "cmo.tarn_r"
        elif(k == 3):
            lvl = np.linspace(-10000, 10000, 127)
            cm = "cmo.tarn_r"

        
        if(usePcolor):
            cp = plt.pcolormesh(
                np.squeeze(x),
                np.squeeze(y),
                np.squeeze(dataPlot),
                cmap=cm,
                vmin=np.min(lvl),
                vmax=np.max(lvl),
            )
        else:
            cp = plt.contourf(
                np.squeeze(x),
                np.squeeze(y),
                np.squeeze(dataPlot),
                lvl,
                extend="both",
                cmap=cm,
            )
        cbar = plt.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        plt.xlabel('Along Fjord [m] %.3f %.3f ∆ Sum: %i' %(np.nanmin(dataPlot),np.nanmax(dataPlot),np.sum(dataPlot)))
        plt.ylabel('Across Fjord [m]')
        plt.title("%s at %.02f days" % (name[k],i/86400.0*dt1))
        if(showZeros):
            cc = plt.contour(
                np.squeeze(x),
                np.squeeze(y),
                np.squeeze(dataPlot),
                [0],
                colors='gray',
                linewidths=0.5
            )
            plt.clabel(cc, inline=3, fontsize=8)
        if(showQuiver):
            u = np.squeeze(dataQuiv[2, zSlice, :, :])
            v = np.squeeze(dataQuiv[4, zSlice, :, :])
            plt.quiver(
                np.squeeze(x),
                np.squeeze(y),
                u/np.sqrt(u**2 + v**2 + 1e-12),
                v/np.sqrt(u**2 + v**2 + 1e-12),
                alpha=.5
                )
        if(isBerg1 or isBerg2):
            cp2 = plt.contourf(np.squeeze(x),
                np.squeeze(y),
                np.squeeze(openFrac[zSlice, :, :]),
                [.1,.5,.9],
                extend="min",
                alpha=.1,
                cmap='cmo.gray')
            #cbar2 = plt.colorbar(cp2)
            #cbar2.set_label('Ocean Fraction')

        j = i/sizeStep + startStep
        str = "figs/compareMelt%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png', dpi=plotDPI)
        plt.close()
        #plt.show()

    os.system('magick -delay %f figs/compareMelt%s*.png -colors 256 -depth 256 figs/compareMelt%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

    #Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/compareMelt*.png')