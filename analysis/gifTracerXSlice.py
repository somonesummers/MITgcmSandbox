from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
import argparse

parser = argparse.ArgumentParser(description='Plot dynamics at xCrossSection')
parser.add_argument('xCrossSection', nargs='?', const=0.0, type=float,
                    help='optional x location [m]')
args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 100
cleanPNGs = True
showDensity = True
usePcolor = False

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
print('Plot DPI:',plotDPI,'; clean PNGs?',cleanPNGs)

if(args.xCrossSection != None):
    print('** Manual xCrossSection detected **')
    xCrossSection = args.xCrossSection

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

maxStep = 0
sizeStep = 1e10
startStep = 1e10

#Decide if iceBerg data files exist
if(os.path.isfile('input/bergMask.bin')):
    isBerg = True
    print('Found icebergs for this run')
else:
    isBerg = False

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

os.system('rm -f figs/sideX*.png')
# os.system('rm -f figs/autosideX*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))

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
    maxDepth = np.zeros(np.shape(x))

    ## deepest contour
    # for i in range(np.shape(x)[1]):
    #     for j in range(np.shape(x)[0]):
    #         bergCount = int(bergsPerCell[j,i])
    #         if(bergMask[j,i] == 1 and bergCount > 0):  #only go in if bergs here
    #             depthFile = 'input/iceberg_depth_%05i.txt' % int(bergMaskNums[j,i])
    #             depths = np.zeros(bergCount)
    #             with open(depthFile,'r') as readFile:
    #                 ii = 0
    #                 for line in readFile:
    #                     if ii >= bergCount:
    #                         print('berg count mismatch in depth')
    #                         break
    #                     depths[ii] = float(line)
    #                     ii += 1
    #             readFile.close()
    #             maxDepth[j,i] = np.max(depths)


    # contourf plot
    openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1


xSlice = np.argmin(np.abs(x[0,:] - xCrossSection))
print('cross section is x =', x[0,xSlice],'index', xSlice)


dynName = ['ptraceDiag','ptraceDiag','ptraceDiag','ptraceDiag']
name = ['TracerTotal','TracerFrac','PlumeAdvect','BergAdvect']
cbarLabel = ["[Vol Frac]","[- Bergs/ + Plume]",'[tracer-m/s]','[tracer-m/s]']

for k in range(len(name)):
    print('\t',name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        lvl = np.linspace(1,1,101)
        cm = "cmo.balance"
        tracerConc = data[0,:,:,:] + data[1,:,:,:] + 1e-8
        tracerFrac = (((data[0,:,:,:]) / tracerConc) * 2) - 1 # -1 all Plume, +1 all Melt
        tracerFrac[tracerConc < 1e-7] = 0
        tracerConc[tracerConc < 1e-7] = 0
        dataVel = mds.rdmds("results/dynDiag", i)
        if k == 0:
            plotData = tracerConc
            lvl = np.linspace(0,.1,101)
            cm = "cmo.turbid"
        elif k == 1: 
            plotData = tracerFrac
            lvl = np.linspace(-1,1,101)
            cm = "cmo.diff"
        elif k == 2: 
            #["Temp", "Sal", "U", "W", "V"]
            plotData = data[0,:,:,:] * dataVel[2,:,:,:]
            lvl = np.linspace(-.005,.005,101)
            cm = "cmo.balance"
        elif k == 3: 
            #["Temp", "Sal", "U", "W", "V"]
            plotData = data[1,:,:,:] * dataVel[2,:,:,:]
            lvl = np.linspace(-.005,.005,101)
            cm = "cmo.balance"
        if(usePcolor):
            cp = plt.pcolormesh(
                np.squeeze(y[:,xSlice]),
                np.squeeze(z),
                np.squeeze(plotData[ :, :, xSlice]),
                cmap=cm,
                vmin=np.min(lvl),
                vmax=np.max(lvl),
            )
        else:
            cp = plt.contourf(
                np.squeeze(y[:,xSlice]),
                np.squeeze(z),
                np.squeeze(plotData[:, :, xSlice]),
                lvl,
                extend="both",
                cmap=cm,
            )
        if(k == 1 or k == 2):
            cc0 = plt.contour(
                    np.squeeze(y[:,xSlice]),
                    np.squeeze(z),
                    np.squeeze(dataVel[2,:, :, xSlice]),
                    [0],
                    colors='white',
                    linestyles='--',
                    linewidths=2,
                    alpha=0.5
                )
        cc = plt.contour(
                    np.squeeze(y[:,xSlice]),
                    np.squeeze(z),
                    np.squeeze(plotData[:, :, xSlice]),
                    [0],
                    colors='gray',
                    linewidths=0.5,
                    alpha=0.5
                )
        plt.clabel(cc, inline=3, fontsize=8)
        plt.plot(y[:,xSlice],topo[:,xSlice],color='black')
        # plt.plot(y[:,xSlice],ice[:,xSlice],color='gray')
        cbar = plt.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        plt.xlabel('Across Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(plotData[ :, :, xSlice]),np.nanmax(plotData[ :, :, xSlice]),np.max(np.isnan(plotData[ :, :, xSlice]))))
        plt.ylabel('Depth [m]')
        plt.title("%s x = %i at %.02f days" % (name[k], x[0,xSlice], i/86400.0*dt))
        j = i/sizeStep + startStep
        
        str = "figs/sideX%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png', dpi=plotDPI)
        # plt.show()
        plt.close()

    if(args.xCrossSection != None):
        os.system('magick -delay %f figs/sideX%s*.png -colors 256 -depth 256 figs/autosideX_%i%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], args.xCrossSection, name[k]))
    else:
        os.system('magick -delay %f figs/sideX%s*.png -colors 256 -depth 256 figs/autosideX_%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/sideX*.png')
