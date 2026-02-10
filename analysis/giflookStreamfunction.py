from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
import argparse

parser = argparse.ArgumentParser(description='Plot streamfunctions width averaged')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='quick option for last from only, double to show plot(s)')
parser.add_argument('-k','--kValues', nargs='*', type=int, default = [0],
                    help='option specification of views to plot [default = stream only]')
parser.add_argument('-n','--numFrames', nargs='?', type=int, default = 60,
                    help='optional specification of numFrames [default = 60]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('-s','--shadow', action='count', default=0,
                    help='option of shadow for mélange [default (off)]')
parser.add_argument('-xl','--xlimit', nargs='?', type=int, default = None,
                    help='optional max x [default = entire frame]')
parser.add_argument('-dpi','--dpi', nargs='?', type=int, default = 125,
                    help='dpi spec [default = 125]')
args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
cleanPNGs = True
showQuiver = True
showZeros = True
showDensity = False
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

print('Plot DPI:',args.dpi,'; clean PNGs:',cleanPNGs, '; usePcolor:', usePcolor)

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

## Custom timeshots available here
# startStep = int(((maxStep - startStep)/startStep) * 1 / 2 )* startStep 
# maxStep = maxStep - 10 * sizeStep

if(args.timeRange != None):
    startStep = args.timeRange[0] * 86400 / dt
    maxStep = args.timeRange[1] * 86400 / dt

if((maxStep-startStep)/sizeStep > args.numFrames):   #if more than numFrames, downscale to be less than numFrames
    dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/args.numFrames)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

#Decide if iceBerg data files exist
if(os.path.isfile(f'results/BRGFlx.{int(startStep):010d}.data') or 
    os.path.isfile(f'results/BRGFlx.{int(startStep):010d}.001.001.data')):
    isBerg = True
    print('Found icebergs for this run')
else:
    isBerg = False
    print('NO icebergs')
# os.system('rm -f figs/side_*.png')
# os.system('rm -f figs/autoside_*.gif')
# os.system('rm -f figs/autoside_*.mov')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))
dz = np.load("input/dz.npy") # we now save this upon grid generation, much better
dy = y[0,0] * 2

if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))

oldAveraging = True
if(isBerg):
    openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1
                if topo[j,i] == 0 and i > 1:
                    openFrac[:,j,i] = 0 #zero weight non-ocean cell
    # print(np.nansum(openFrac[:,:,0:3],axis=(1)))
else:
    openFrac = mds.rdmds("results/hFacC")
    openFrac[:,:,0] = 1
# print(np.sum(openFrac[:,:,:],axis=(0,1)))
if(np.min(openFrac) < 0 or np.max(openFrac) > 1):
    print(f'Error in openFrac min/max: {np.min(openFrac):0.2f}/{np.max(openFrac):0.2f}')
print('averaging over all cross sections')


dynName = ['dynDiag','dynDiag','dynDiag']
name = ["Psi", "UfromPsi", "WfromPsi"]
cbarLabel = ["[mSv]","[m/s]","[m/s]"]


if(args.quick > 0):
    startStep = maxStep
    cleanPNGs = False
if(args.kValues == None):
    kList = range(len(name))
else:
    kList = args.kValues
for k in kList:
    print("\t" + name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        if(showQuiver):
            dataQuiv = mds.rdmds("results/dynDiag", i)
        if(isBerg):
            dataBergs = mds.rdmds("results/BRGFlx",i)
            if(dataBergs.shape[0] == 6):
                openFrac[:,:,:] = dataBergs[5,:,:,:] #directly saved for these runs
                openFrac[:,topo==0] = 0
                openFrac[:,:,0] = 1

                oldAveraging = False
        if(args.shadow > 0): #enable berg shadows here
            if(dataBergs.shape[0] < 6):
                dataBergPlot = dataBergs[0,:,:,:]
                dataBergPlot[dataBergPlot != 0] = .2
            else:
                dataBergPlot = dataBergs[5,:,:,:]
            dataBergPlot[dataBergPlot == 0] = 1 #np.nan
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        kk = k
        if k == 0:
            lvl = np.linspace(-45,45,37)
            cm = "bwr"
            kk = 2
            dataPsi = -1 * np.cumsum(np.average(data[kk, ::-1, :, :],axis=1,weights=openFrac[::-1,:,:]) * dy * np.sum(openFrac[::-1,:,:],axis=1) * dz[::-1,None],axis=0)[::-1,:] * 1e-3
            # dataPsi = np.cumsum(np.sum(data[kk, :, :, :],axis=1) * (dy),axis=0) * 1e-3
        if k == 1:
            lvl = uRange
            cm = uCmap
            kk = 2
            dataPsi = -1 * np.cumsum(np.average(data[kk, ::-1, :, :],axis=1,weights=openFrac[::-1,:,:]) * dy * np.sum(openFrac[::-1,:,:],axis=1) * dz[::-1,None],axis=0)[::-1,:] * 1e-3
            dataPsi = 1 * np.gradient(dataPsi,axis=0) / (dy * np.sum(openFrac,axis=1)) / dz[:,None] * 1e3
            # dataPsi = np.cumsum(np.sum(data[kk, :, :, :],axis=1) * (dy),axis=0) * 1e-3
        if k == 2:
            lvl = wRange
            cm = wCmap
            kk = 2
            dataPsi = -1 * np.cumsum(np.average(data[kk, ::-1, :, :],axis=1,weights=openFrac[::-1,:,:]) * dy * np.sum(openFrac[::-1,:,:],axis=1) * dz[::-1,None],axis=0)[::-1,:] * 1e-3
            dataPsi = 1 * np.gradient(dataPsi,axis=1) / (dy * np.sum(openFrac,axis=1)) / dz[:,None] * 1e3
            # dataPsi = np.cumsum(np.sum(data[kk, :, :, :],axis=1) * (dy),axis=0) * 1e-3
        plt.figure(figsize=(12, 5))
        if(usePcolor):
            cp = plt.pcolormesh(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(dataPsi),
                cmap=cm,
                vmin=np.min(lvl),
                vmax=np.max(lvl),
            )
        else:
            cp = plt.contourf(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(dataPsi),
                lvl,
                extend="both",
                cmap=cm,
            )
        # if(k == 5): #special bounds to highlight refreezing areas
        #     cp.cmap.set_under('r')
        plt.plot(x[0,:],topo[int(np.shape(x)[0]/2),:],color='black')
        if(args.shadow > 0):
            cp2 = plt.contourf(
                x[0,:],
                np.squeeze(z),
                np.squeeze(np.nanmean(dataBergPlot[:, 1:-1, :],axis=1)),
                [.4,.6,.8,.9,.95],
                extend="min",
                alpha=.2,
                cmap='cmo.gray')
            #cbar2 = plt.colorbar(cp2)
            #cbar2.set_label('Ocean Fraction')
        cbar = plt.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        if(k > 0 and (dynName[k] == 'dynDiag')):
            salt = np.squeeze(np.average(data[1,:,1:-1,:],weights=openFrac[:,1:-1,:],axis=1))
            if(i == startStep): #only calc pressure once
                pressure = -1 * np.ones(salt.shape) * 1020 * 9.81 * np.repeat(np.expand_dims(z,1), salt.shape[1], axis=1) /10e3
            CT = gsw.CT_from_t(salt, np.average(data[0,:,1:-1,:],weights=openFrac[:,1:-1,:],axis=1), pressure)
            density = gsw.rho(salt, CT, 0) - 1000 #in-stu density less 1000
            densityLevels = np.linspace(25,28,16)
            cc = plt.contour(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(density),
                densityLevels,
                colors='black',
                linewidths=0.5,
                alpha=0.5
            )
            plt.clabel(cc, inline=3, fontsize=8)
            c0 = plt.contour(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(dataPsi),
                [0],
                colors='gray',
                linewidths=0.5,
                alpha=.5
                )
            plt.clabel(c0, inline=3, fontsize=8)
        if(showZeros and k == 0):
            cc = plt.contour(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(dataPsi),
                lvl[2:-2:4],
                colors='black',
                linewidths=0.5,
                alpha=1
            )
            plt.clabel(cc, inline=3, fontsize=8)
        if(showQuiver):
            downSampleX = 5
            downSampleZ = 1
            # np.average(data[kk, ::-1, :, :],axis=1,weights=openFrac[::-1,:,:])
            u = np.squeeze(np.average(dataQuiv[2, :, :, :],axis=1,weights = openFrac))[::downSampleZ,::downSampleX]
            w = np.squeeze(np.average(dataQuiv[3, :, :, :],axis=1,weights = openFrac))[::downSampleZ,::downSampleX]
            plt.quiver(
                np.squeeze(x)[0,::downSampleX],
                np.squeeze(z)[::downSampleZ],
                u,#/np.sqrt(u**2 + w**2 + 1e-12),
                w,#/np.sqrt(u**2 + w**2 + 1e-12),
                alpha=.2,
                width = .002, #width of line
                scale = .5,
                angles='xy',#60
                pivot='mid'
                )
        if(args.xlimit != None):
            plt.xlim([0, args.xlimit])
        plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(dataPsi),np.nanmax(dataPsi),np.max(np.isnan(dataPsi))))
        plt.ylabel('Depth [m]')
        plt.title("%s at %.02f days" % (name[k], i/86400.0*dt))
        j = i/sizeStep
        
        str = "figs/sideStream_%s%05i.png" % (name[k],j)
        plt.tight_layout()
        plt.savefig(str, format='png', dpi=args.dpi)
        if(args.quick > 1):
            plt.show()
        plt.close()

    
    if(args.quick == 0):
        timeStr = ''
        if(args.timeRange != None):
            timeStr = f'_T_{args.timeRange[0]}_{args.timeRange[1]}'
        os.system('magick -delay %f figs/sideStream_%s*.png -depth 32 +dither -verbose figs/sideStream_%s%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k],timeStr))
        if(makeMovie):
            os.system('ffmpeg -r %f -i figs/sideStream_%s%%05d.png -c:v libx264 -r 30 -pix_fmt yuv420p figs/sideStream_%s.mov' %(80/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/sideStream_*.png')