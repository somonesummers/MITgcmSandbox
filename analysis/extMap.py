from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
import argparse

parser = argparse.ArgumentParser(description='Plot dynamics at zDepth')
parser.add_argument('-z','--zDepth', nargs=1, type=float,default = None,
                    help='optional depth location [m]')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='quick option for last from only, double to show plot(s)')
parser.add_argument('-k','--kValues', nargs='*', type=int, default = None,
                    help='option specification of views to plot [default = all]')
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
parser.add_argument('--quiver', nargs="?", type=int, default = 1,
                    help='specify if quiver arrows are shown [Default = True]')
args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
cleanPNGs = True
showQuiver = False
showZeros = True
showDensity = True
makeMovie = False   

showQuiver = args.quiver
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

if(args.zDepth != None):
    print('** Manual zDepth detected **')
    zDepth = args.zDepth

os.system('rm -f figs/map*.png')

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
if(os.path.isfile('input/bergMask.bin')):
    isBerg = True
    print('Found icebergs for this run')
else:
    isBerg = False

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))

if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))

oldAveraging = True
if(isBerg):
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    bergMaskNums = np.fromfile('input/bergMaskNums.bin', dtype='>f8')
    bergMaskNums = bergMaskNums.reshape(np.shape(x))
    bergsPerCell = np.fromfile('input/numBergsPerCell.bin', dtype='>f8')
    bergsPerCell = bergsPerCell.reshape(np.shape(x))
    bergContaingingCells = int(np.sum(bergMask))
    maxDepth = np.zeros(np.shape(x))

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
                    
if(np.min(openFrac) < 0 or np.max(openFrac) > 1):
    print(f'Error in openFrac min/max: {np.min(openFrac):0.2f}/{np.max(openFrac):0.2f}')


zSlice = np.argmin(np.abs(z[:]- zDepth))
print('depth is z =', z[zSlice], 'index', zSlice)

# dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'BRGFlx','ptraceDiag','ptraceDiag']
# name = ["Temp", "Sal", "U", "W", "V", "BRGmltRt",'TracePlume','TraceBerg']
# cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[m/d]","[Vol Frac]","[Vol Frac]"]

dynNameList = ['heatDiag','momDiag','presDiag']
# dynNameList = ['presDiag']
for dynName in dynNameList:
    if(dynName == 'presDiag'):
        name = ['PHIHYD', 'PHI_NH']
    elif(dynName == 'momDiag'):
        name = ['Um_Diss','Um_Advec','Um_Cori','Um_dPhiX','Wm_Diss','Wm_Advec']
    elif(dynName == 'heatDiag'):
        name = ['ADVx_TH','ADVy_TH','ADVr_TH','UTHMASS','VTHMASS','WTHMASS']
    # print(metadata)
    if(args.quick > 0):
        startStep = maxStep
        cleanPNGs = False
    if(args.kValues == None):
        kList = range(len(name))
    else:
        kList = args.kValues
    for k in kList:
        print("\t" + dynName + ": " + name[k])
        for i in np.arange(startStep, maxStep + 1, sizeStep):
            if(showQuiver or showDensity):
                dataQuiv = mds.rdmds("results/dynDiag", i)
            if(isBerg):
                dataBergs = mds.rdmds("results/BRGFlx",i)
                if(dataBergs.shape[0] == 6):
                    # print('ttest')
                    openFrac[:,:,:] = dataBergs[5,:,:,:] #directly saved for these runs
                    for jj in range(np.shape(x)[0]): #clean up non-berg parts of this mask
                        for ii in range(np.shape(x)[1]):
                            if bergMask[jj,ii] == 0:
                                openFrac[:,jj,ii] = 1
                            if topo[jj,ii] == 0 and ii > 1:
                                openFrac[:,jj,ii] = 0 #zero weight non-ocean cell
                    # print(np.sum(openFrac[1,:,:],axis=0))
                    # print(np.sum(openFrac2[1,:,:],axis=0))
                    # oldAveraging = False
            # print(np.sum(openFrac[:,:,:],axis=(0,1)))
            if(args.shadow > 0): #enable berg shadows here
                if(dataBergs.shape[0] < 6):
                    dataBergPlot = dataBergs[0,:,:,:]
                    dataBergPlot[dataBergPlot != 0] = .2
                else:
                    dataBergPlot = dataBergs[5,:,:,:]
                dataBergPlot[dataBergPlot == 0] = 1 #np.nan
            data = mds.rdmds("results/%s"%(dynName), i)
            kk = k
            if(dynName == 'momDiag'):
                cm = 'cmo.balance'
                lvl = np.linspace(-5e-5,5e-5,31)
            elif(name[k] == 'PHI_NH'):
                cm = 'PuOr_r'
                lvl = np.linspace(-.4,.4,31)
            elif(name[k] == 'PHIHYD'):
                cm = 'PuOr_r'
                data[:,dataQuiv[2,:,:,:] == 0] = np.nan
                minVal = np.nanmin(data[k,zSlice,:,:])
                maxVal = np.nanmax(data[k,zSlice,:,:])
                mn = np.nanmean(data[k,zSlice,:,2:-2])
                data[k,:,:,:] -= mn
                sd = np.nanstd(data[k,zSlice,:,:])
                lvl = np.linspace(-.4,.4,31)
            else:
                minVal = np.nanmin(data[k,zSlice,:,:])
                maxVal = np.nanmax(data[k,zSlice,:,:])
                mn = np.nanmean(data[k,zSlice,:,:])
                sd = np.nanstd(data[k,zSlice,:,:])
                if((maxVal)/(minVal) < .1 ):
                    cm = 'cmo.balance'
                    lvl = np.linspace(-mn-3*sd,mn+3*sd,31)
                else:
                    cm = 'cmo.tempo'
                    lvl = np.linspace(mn-1*sd,mn+1*sd,31)
            plt.figure(figsize=(12, 5))
            dataPlot = data[kk, zSlice, :, :]
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
            cbar = plt.colorbar(cp,orientation="horizontal",fraction=0.06,format='%.2f')
            # cbar.set_label(cbarLabel[k])
            ax1 = plt.gca()
            ax1.set_aspect('equal')
            plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(dataPlot),np.nanmax(dataPlot),np.max(np.isnan(dataPlot))))
            plt.ylabel('Across Fjord [m]')
            plt.title("%s depth %.2f m at %.02f days" % (name[k], z[zSlice] ,i/86400.0*dt))
            plt.contour(x,
                        y,
                        topo,
                        [zDepth],
                        colors = 'black')
            if(args.shadow > 0):
                cp2 = plt.contourf(np.squeeze(x),
                np.squeeze(y),
                np.squeeze(dataBergPlot[zSlice, :, :]),
                [.4,.6,.8,.9,.95],
                extend="min",
                alpha=.1,
                cmap='cmo.gray')
                #cbar2 = plt.colorbar(cp2)
                #cbar2.set_label('Ocean Fraction')
            if(showQuiver):
                downSampleX = 5
                downSampleY = 3
                u = np.squeeze(dataQuiv[2, zSlice, :, :])[::downSampleY,::downSampleX]
                v = np.squeeze(dataQuiv[4, zSlice, :, :])[::downSampleY,::downSampleX]
                plt.quiver(
                    np.squeeze(x)[::downSampleY,::downSampleX],
                    np.squeeze(y)[::downSampleY,::downSampleX],
                    u/np.sqrt(u**2 + v**2 + 1e-12),
                    v/np.sqrt(u**2 + v**2 + 1e-12),
                    alpha=.2,
                    width = .001, #width of line
                    scale = 120
                    )
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
            if(args.xlimit != None):
                plt.xlim([0, args.xlimit])
            j = i/sizeStep
            
            str = "figs/Map_%s%05i.png" % (name[k],j)
            plt.tight_layout()
            plt.savefig(str, format='png', dpi=args.dpi)
            if(args.quick > 1):
                plt.show()
            plt.close()

        
        if(args.quick == 0):
            os.system('magick -delay %f figs/Map_%s*.png -colors 256 -depth 256 figs/autoMap_%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))
            if(makeMovie):
                os.system('ffmpeg -r %f -i figs/Map_%s%%05d.png -c:v libx264 -r 30 -pix_fmt yuv420p figs/autoMap_%s.mov' %(80/((maxStep-startStep)/sizeStep), name[k], name[k]))

    #Clean up intermediate pngs
        if(cleanPNGs):
            os.system('rm -f figs/Map_*.png')