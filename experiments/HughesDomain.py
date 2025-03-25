from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
# import argparse

# parser = argparse.ArgumentParser(description='Plot dynamics at ySlice')
# parser.add_argument('yCrossSection', nargs='?', const=0.0, type=float,
#                     help='optional slice location [m]')
# args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 100
cleanPNGs = True
showQuiver = False
showZeros = True
showDensity = True
makeMovie = False   

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

# if(args.yCrossSection != None):
#     print('** Manual ySlice detected **')
#     yCrossSection = args.yCrossSection
# #Overwrite local settings here if desired
# # usePcolor = True

print('Plot DPI:',plotDPI,'; clean PNGs:',cleanPNGs, '; usePcolor:', usePcolor)

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

# os.system('rm -f figs/side_*.png')
# os.system('rm -f figs/autoside_*.gif')
# os.system('rm -f figs/autoside_*.mov')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))
hFacC = mds.rdmds("results/hFacC")

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
                    

print('averaging over all cross sections')


dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'BRGFlx', 'BRGFlx', 'BRGFlx']
name = ["Temp", "Sal", "u", "W", "V", "BRGmltRt", 'BRG_TauX', 'BRG_TauY']
cbarLabel = ["[C]", "[ppt]", "u [m/s]", "[m/s]", "[m/s]", "[m/d]", "[kN/m^2]", "[kN/m^2]"]

# fig, axs = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
plt.figure(figsize=(12, 6))
ax1 = plt.subplot(2,2,1)
ax2 = plt.subplot(2,4,3)
ax3 = plt.subplot(2,4,4)
ax4 = plt.subplot(2,1,2)
axes = [ax1, ax2, ax3, ax4]
labels = ["A","B","C","D"]
# plt.figure(figsize=(12, 8),layout="constrained")
for k in [0]:
    print("\t" + name[k])
    for i in [maxStep]:
        if(showQuiver):
            dataQuiv = mds.rdmds("results/dynDiag", i)
        if(isBerg and os.path.isfile('results/BRGFlx.%010i.001.001.data' % i)):
            localBergs = True
        else:
            localBergs = False
        if((not localBergs) and k==5):
            #fill melt image with 0s if bergs in run but not frame
            data = np.zeros(np.shape(mds.rdmds("results/%s"%(dynName[k-1]), i)))
        else:
            data = mds.rdmds("results/%s"%(dynName[k]), i)
        
        if k == 0:
            lvl = np.linspace(-0.5,3.0,127)
            cm = tempCmap
        elif k == 1:
            lvl = saltRange
            cm = saltCmap
        elif k == 2:
            lvl = np.linspace(-.15,.15,127)
            cm = uCmap
        elif k == 3:
            lvl = wRange
            cm = wCmap
        elif k == 4:
            lvl = vRange
            cm = vCmap
        elif k == 5:
            lvl = meltRange
            cm = meltCmap
        elif k == 6:
            lvl = np.linspace(-10,10,127)
            cm = "cmo.balance"
        elif k == 7:
            lvl = np.linspace(-10,10,127)
            cm = "cmo.balance"
        if(k > 4):
            kk = k - 3
            if(k > 5):
                data = data/1000 #Pa to kPa
        else:
            kk = k
        if(usePcolor):
            cp = ax4.pcolormesh(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(np.average(data[kk, :, 1:-1, :], weights=openFrac[:,1:-1,:],axis=1)),
                cmap=cm,
                vmin=np.min(lvl),
                vmax=np.max(lvl),
            )
        else:
            cp = plt.contourf(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(np.average(data[kk, :, 1:-1, :], weights=openFrac[:,1:-1,:],axis=1)),
                lvl,
                extend="both",
                cmap=cm,
            )
        plt.plot(x[0,:],topo[1,:],color='black')
        if(localBergs):
            cp2 = plt.contourf(
                x[0,:],
                np.squeeze(z),
                np.squeeze(np.average(1-openFrac[:, 1:-1, :],axis=1)),
                [.05,.1,.2,.4,.8],
                extend="max",
                alpha=.5,
                cmap='cmo.gray_r')
            #cbar2 = plt.colorbar(cp2)
            #cbar2.set_label('Ocean Fraction')
        cbar = plt.colorbar(cp,ticks=np.linspace(lvl[0],lvl[-1],7))
        cbar.set_label(cbarLabel[k])
        if(showZeros):
            if( k == 2 or k == 6 or k == 7):
                cc = plt.contour(
                    np.squeeze(x[0,:]),
                    np.squeeze(z),
                    np.squeeze(np.average(data[kk, :, 1:-1, :], weights=openFrac[:,1:-1,:],axis=1)),
                    [0],
                    colors='gray',
                    linewidths=0.5,
                    alpha=0.5
                )
                plt.clabel(cc, inline=3, fontsize=8)
        cbar = plt.colorbar(cp2,ticks=[.8,.4,.2,.1,0.05])
        cbar.set_label('Fraction Ice')       
        # ax4.set_title("Width Averaged %s for $\\lambda$ = 0.20, U = 0.12 m/s" % (name[k]))
        ax4.set_title("Width Averaged %s after %s Days" %(name[k],maxStep*dt/86400))
        ax4.set_ylabel('Depth [m]')
        ax4.set_xlabel('Along Fjord [m]')
        j = i/sizeStep
        
        # str = "figs/HughesDomain_%s%05i.png" % (name[k],j)
        
        # plt.savefig(str, format='png', dpi=plotDPI)
        # plt.show()
        # plt.close()

# Iceberg depths
bergDepths = np.fromfile('input/icebergs_depths.bin', dtype='>f8')
# bergDepths = bergDepths.reshape((500,14,160))
bergDepths[bergDepths == 0] = np.nan

ax3.hist(bergDepths,bins = 50)
ax3.set_ylabel('Count')
ax3.set_xlabel('Depth [m]')
ax3.set_title('Iceberg Depths')

# Iceberg Coverage
pltHelper = 1-hFacC[0,:,:]
pltHelper[bergMask == 0] = np.nan
pc = ax1.pcolormesh(np.squeeze(x[0,:]),np.squeeze(y[:,0]),pltHelper,cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
ax1.set_title('Iceberg Surface Areal Fraction $\\lambda$')
ax1.set_ylabel('Across Fjord [m]')
ax1.set_xlabel('Along Fjord [m]')
# ax1.set_xlim([10000,18500])
ax1.set_xlim([0,42000])
cbar.set_label('$\\lambda$')

#Icebergs with Depth
hFacC[:,:,hFacC[0,1,:] == 1] = np.nan #nan an no-berg cells
for i in np.where(bergMask[1,:]==1)[0]:
    for j in np.where(bergMask[:,i]==1)[0]:
        ax2.plot(1-hFacC[:,j,i],z,alpha=.5,color='xkcd:gray',linewidth=.5)
ax2.plot(1-np.nanmean(openFrac[:,1:-1,bergMask[1,:]==1],axis=(1,2)),z,alpha=1,color='xkcd:black',linewidth=1,linestyle='--',label='Mean')
hlper = 1-np.nanmean(openFrac[0,1:-1,bergMask[1,:]==1])
ax2.plot([hlper,hlper],[z[0],z[-1]],alpha=.25,color='xkcd:black',linewidth=1,linestyle='--',label='Surface Mean')
ax2.set_title('$\\lambda(z)$')
ax2.set_ylabel('Depth [m]')
ax2.set_xlabel('$\\lambda$')
ax2.legend()


plt.tight_layout()

for label,ax in zip(labels,axes):
    ax.text(
        ax.get_xlim()[0], ax.get_ylim()[1], label,
        fontsize='x-large', va='bottom',ha='right', fontfamily='sans serif')

str = "figs/HughesDomain.png"
plt.savefig(str, format='png', dpi=plotDPI)
plt.show()
plt.close()



        