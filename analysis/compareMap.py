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
if(os.path.isfile('input/bergMask.bin')):
    isBerg = True
    print('Found icebergs for this run')
else:
    isBerg = False

#Clean up old gifs and pngs
os.system('rm -f figs/compareMap*.png')
os.system('rm -f figs/compareMap*.gif')

y = mds.rdmds("results/YC")
x = mds.rdmds("results/XC")
z = mds.rdmds("results/RC")



if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))

if(isBerg):
    # contourf plot
    openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1


zSlice = np.argmin(np.abs(z[:,0,0]- zDepth))
print('depth is z =', z[zSlice,0,0], 'index', zSlice)

# iceEdge = np.interp(z[zSlice,0,0],ice[0,:],x[0,:])

name = ["Temp", "Sal", "U", "W", "V"]
cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]"]

if(isBerg):
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag', 'BRGFlx']
    name = ["Temp", "Sal", "U", "W", "V", "SPD", "BRGmltRt"]
    cbarLabel = ["[∆ C]", "[∆ ppt]", "[∆ m/s]", "[∆ m/s]", "[∆ m/s]","[∆ m/s]", "[∆ m/d]"]
else:
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag']
    name = ["Temp", "Sal", "U", "W", "V", "SPD"]
    cbarLabel = ["[∆ C]", "[∆ ppt]", "[∆ m/s]", "[∆ m/s]", "[∆ m/s]","[∆ m/s]"]

for k in range(len(name)):
    print("\t" + name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        if(showQuiver):
            dataQuiv = mds.rdmds(folder1 + "/dynDiag", i)
        if(isBerg and os.path.isfile(folder1 + '/BRGFlx.%010i.001.001.data' % i)):
            localBergs = True
        else:
            localBergs = False
        if((not localBergs) and dynName[k]=='BRGFlx'):
            #fill melt image with 0s if bergs in run but not frame
            data1 = np.zeros(np.shape(mds.rdmds(folder1 + "/%s"%(dynName[k-1]), i)))
            data2 = np.zeros(np.shape(mds.rdmds(folder2 + "/%s"%(dynName[k-1]), i * dt1 / dt2)))
        else:
            data1 = mds.rdmds(folder1 + "/%s"%(dynName[k]), i)
            data2 = mds.rdmds(folder2 + "/%s"%(dynName[k]), i * dt1 / dt2)
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
            lvl = np.linspace(-0.001, 0.001, 127)
            cm = "cmo.curl"
        elif k == 6: #meltRate
            lvl = np.linspace(-.1, .1, 127)
            cm = "cmo.curl"
        elif k == 5: #speed
            lvl = np.linspace(-.1, .1, 127)
            cm = "cmo.curl"
        if(k == 6):
            kk = 2
        else:
            kk = k

        if(name[k] == "SPD"):
            spd1 = np.sqrt(data1[2, zSlice, :, :]**2 + data1[3, zSlice, :, :]**2 + data1[4, zSlice, :, :]**2)
            spd2 = np.sqrt(data2[2, zSlice, :, :]**2 + data2[3, zSlice, :, :]**2 + data2[4, zSlice, :, :]**2)
            dataPlot = spd1 - spd2
        else:
            dataPlot = data1[kk, zSlice, :, :] - data2[kk, zSlice, :, :]
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
        plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(dataPlot),np.nanmax(dataPlot),np.max(np.isnan(dataPlot))))
        plt.ylabel('Across Fjord [m]')
        plt.title("%s depth %f at %.02f days" % (name[k], z[zSlice,0,0] ,i/86400.0*dt1))
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
        if(localBergs):
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
        str = "figs/compareMap%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png', dpi=plotDPI)
        plt.close()
        #plt.show()

    os.system('magick -delay %f figs/compareMap%s*.png -colors 256 -depth 256 figs/compareMap%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

    #Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/compareMap*.png')