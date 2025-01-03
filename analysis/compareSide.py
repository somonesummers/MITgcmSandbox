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
zDepth = 50
plotDPI = 100
cleanPNGs = True
showZeros = True


if(os.path.isfile('input/plotHelperLocal.py')):
    sys.path.append('input')
    from plotHelperLocal import *
    print('Found experiment plotting settings')
elif(os.path.isfile('../plotHelper.py')):
    print('no custom plotting settings, using local default')
    from plotHelper import *
else:  
    print('no defaults found')

print('Plot DPI:',plotDPI,'; clean PNGs:',cleanPNGs, '; usePcolor:', usePcolor)

maxStep = 0
sizeStep = 1e10
startStep = 1e10

for file in os.listdir(folder1):
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

if(maxStep/sizeStep > 50):   #if more than 50 frames, downscale to be less than 50
    dwnScale = round((maxStep/sizeStep)/50)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)
#Decide if iceBerg data files exist
if(os.path.isfile('input/data.iceberg')):
    isBerg = True
else:
    isBerg = False

os.system('rm -f figs/compareSide_*.png')
os.system('rm -f figs/autoCompareSide_*.gif')

interpolateGrid = False
x = mds.rdmds("%s/XC" % folder1)
y = mds.rdmds("%s/YC" % folder1)
z = mds.rdmds("%s/RC" % folder1)

x2 = mds.rdmds("%s/XC" % folder2)
y2 = mds.rdmds("%s/YC" % folder2)
z2 = np.squeeze(mds.rdmds("%s/RC" % folder2))

if(np.max(x != x2)):
    print("\tUneven lateral grid detected, attempting to interpolate onto %s grid" %folder1)
    interpolateGrid = True
if(np.max(z != z2)):
    print("\tUneven vertical grid detected, attempting to interpolate onto %s grid" %folder1)
    interpolateGrid = True

if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
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

    # contourf plot
    openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1
                    

ySlice = np.argmin(np.abs(y[:,0] - yCrossSection))
print('cross section is y =', y[ySlice,0], 'index', ySlice)

# if(isBerg):
#     dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'BRGFlx']
#     name = ["Temp", "Sal", "U", "W", "V","BRGmltRt"]
#     cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[m/d]"]
# else:
dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag','dynDiag','BRGFlx']
name = ["Temp", "Sal", "U", "W", "V","SPD_h","BRGmltRt"]
cbarLabel = ["[∆ C]", "[∆ ppt]", "[∆ m/s]", "[∆ m/s]", "[∆ m/s]","[∆ m/s]","[∆ m/d]"]

for k in range(len(name)):
    print('\t',name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data1 = mds.rdmds("%s/%s" % (folder1, dynName[k]), i)
        data2 = mds.rdmds("%s/%s" % (folder2, dynName[k]), i * dt1/dt2)
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
        elif k == 5: #speed
            lvl = np.linspace(-.1, .1, 127)
            cm = "cmo.curl"
        elif k == 6: #meltRate
            lvl = np.linspace(-.1, .1, 127)
            cm = "cmo.curl"
        if(k == 6):
            kk = 2
        else:
            kk = k
        if(k == 5): #speed plot
            if(interpolateGrid):
                from scipy.interpolate import RegularGridInterpolator as rgi
                data_interp_u = rgi((z2,y2[:,0],x2[0,:]), data2[2,:,:,:],bounds_error=False)
                data_interp_v = rgi((z2,y2[:,0],x2[0,:]), data2[4,:,:,:],bounds_error=False)
                points = np.array(np.meshgrid(z,y[:,0],x[0,:])).T.reshape(-1, 3)
                data_new_u = data_interp_u(points)
                data_new_v = data_interp_v(points)
                shp = np.shape(data1[2,:,:,:])
                data_new_u = data_new_u.reshape([shp[2],shp[0],shp[1]]).T
                data_new_u = np.swapaxes(data_new_u,0,1) #kind of a mess, but this gets it back into z,y,x order
                data_new_v = data_new_v.reshape([shp[2],shp[0],shp[1]]).T
                data_new_v = np.swapaxes(data_new_v,0,1) #kind of a mess, but this gets it back into z,y,x order
                data = np.sqrt(data1[2,:,:,:]**2 + data1[4,:,:,:]**2) - np.sqrt(data_new_u[:,:,:]**2 + data_new_v[:,:,:]**2)
            else:
                data = np.sqrt(data1[2,:,:,:]**2 + data1[4,:,:,:]**2) - np.sqrt(data2[2,:,:,:]**2 + data2[4,:,:,:]**2)
        else:
            if(interpolateGrid):
                from scipy.interpolate import RegularGridInterpolator as rgi
                data_interp = rgi((z2,y2[:,0],x2[0,:]), data2[kk,:,:,:],bounds_error=False)
                points = np.array(np.meshgrid(z,y[:,0],x[0,:])).T.reshape(-1, 3)
                data_new = data_interp(points)
                shp = np.shape(data1[kk,:,:,:])
                data_new = data_new.reshape([shp[2],shp[0],shp[1]]).T
                data_new = np.swapaxes(data_new,0,1) #kind of a mess, but this gets it back into z,y,x order
                data = data1[kk,:,:,:] - data_new[:,:,:]
            else:
                data = data1[kk,:,:,:] - data2[kk,:,:,:]
        if(usePcolor):
            cp = plt.pcolormesh(
                np.squeeze(x[ySlice,:]),
                np.squeeze(z),
                np.squeeze(data[ :, ySlice, :]),
                cmap=cm,
                vmin=np.min(lvl),
                vmax=np.max(lvl),
            )
        else:
            cp = plt.contourf(
                np.squeeze(x[ySlice,:]),
                np.squeeze(z),
                np.squeeze(data[ :, ySlice, :]),
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
            # cbar2 = plt.colorbar(cp2)
            # cbar2.set_label('Ocean Fraction')
        if(showZeros):
            cc = plt.contour(
                np.squeeze(x[ySlice,:]),
                np.squeeze(z),
                np.squeeze(data[ :, ySlice, :]),
                [0],
                colors='gray',
                linewidths=0.5
            )
            plt.clabel(cc, inline=3, fontsize=8)
        cbar = plt.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(data[ :, ySlice, :]),np.nanmax(data[ :, ySlice, :]),np.max(np.isnan(data[ :, ySlice, :]))))
        plt.ylabel('Depth [m]')
        plt.title("%s y = %i at %i days" % (name[k], y[ySlice,0], i/86400*dt1))
        j = i/startStep
        
        str = "figs/compareSide_%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png')
        plt.close()
        #plt.show()

    os.system('magick -delay %f figs/compareSide_%s*.png -colors 256 -depth 256 figs/autoCompareSide_%s.gif' %(500/(maxStep/sizeStep), name[k], name[k]))
    if(cleanPNGs):
        os.system('rm -f figs/compareSide_*.png')
