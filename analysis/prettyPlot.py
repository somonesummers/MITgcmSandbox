from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
import argparse
import glob
import pickle
from scipy.integrate import simpson
sys.path.append('/Users/psummers8/Documents/glaciome1D')
from glaciome1D import constants, glaciome

parser = argparse.ArgumentParser(description='Plot multifigure pretty plot')
parser.add_argument('-i','--iteration', nargs='?', type=int, default = None,
                    help='optional specification iteration to plot [default = None]')
parser.add_argument('-f','--files', nargs='?', default='couplingResults/MITgcmRun_',
                    help='file string used to name pickle files FILES00000.pickle [default = couplingResults/MITgcmRun_]')
parser.add_argument('--noDensity', action='count', default=0,
                    help='quick option to disable density levels')
parser.add_argument('--noZero', action='count', default=0,
                    help='quick option to disable 0 contour in speed')
parser.add_argument('--noBergs', action='count', default=0,
                    help='quick option to disable berg contours')
args = parser.parse_args()


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
usePcolor = False
showDensity = True
showZeros = True
showBergs = True
plotDPI = 300

if(args.noDensity > 0):
    showDensity = False
if(args.noZero > 0):
    showZeros = False
if(args.noBergs > 0):
    showBergs = False

print('Plot DPI:',plotDPI,'; usePcolor:', usePcolor)

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
sizeStep = sizeStep*10

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

mitgcmTime = np.zeros(int(maxStep/sizeStep))
mitgcmTemp = np.zeros(int(maxStep/sizeStep))
for i in range(int(maxStep/sizeStep)):
    mitgcmTime[i] = (i+1)*dt*sizeStep
    data = mds.rdmds("results/%s"%('dynDiag'), (i+1)*sizeStep)
    mitgcmTemp[i] = np.max(data[0,:,:,:])


#Decide if iceBerg data files exist
if(os.path.isfile('input/bergMask.bin')):
    isBerg = True
    print('Found icebergs for this run')
else:
    isBerg = False


# os.system('rm -f figs/side_*.png')
# os.system('rm -f figs/autoside_*.gif')
# os.system('rm -f figs/autoside_*.mov')

x = mds.rdmds("results/XC")/1e3
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
                    

dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'BRGFlx','ptraceDiag','ptraceDiag']
name = ["Temperature", "Sal", "Along Fjord Velocity U", "W", "V", "BRGmltRt",'TracePlume','TraceBerg']
cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[m/d]","[Vol Frac]","[Vol Frac]"]

fig, axes = plt.subplots(2, 2, figsize=(14, 6), layout="constrained",width_ratios=[3, 1])
ax1 = axes[0,0]
ax2 = axes[1,0]
ax3 = axes[0,1]
ax4 = axes[1,1]

if(args.iteration == None):
    i = maxStep
else:
    i = args.iteration

for ax, k in zip([ax1,ax2],[0,2]):
    data = mds.rdmds("results/%s"%(dynName[k]), i)
    dataBerg = mds.rdmds("results/BRGFlx", i)
    kk = k
    if k == 0:
        lvl = np.linspace(-2,5.5,21)
        cm = tempCmap
    elif k == 1:
        lvl = saltRange
        cm = saltCmap
    elif k == 2:
        lvl = np.linspace(-.2,.2,21)
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
        kk = k - 3
    elif k == 6:
        lvl = np.linspace(0,0.05,127)
        cm = "cmo.matter"
        kk = 0
    elif k == 7:
        lvl = np.linspace(0,0.05,128)
        cm = "cmo.matter"
        kk = 1

    plotData = np.average(data[kk, :, 1:-1, :], weights=openFrac[:,1:-1,:],axis=1)
    if(usePcolor):
        cp = ax.pcolormesh(
            np.squeeze(x[0,:]),
            np.squeeze(z),
            np.squeeze(plotData),
            cmap=cm,
            vmin=np.min(lvl),
            vmax=np.max(lvl),
        )
    else:
        cp = ax.contourf(
            np.squeeze(x[0,:]),
            np.squeeze(z),
            np.squeeze(plotData),
            # np.squeeze(data[kk, :, 7, :]),
            lvl,
            extend="both",
            cmap=cm,
        )
    ax.plot(x[0,:],topo[int(np.shape(x)[0]/2),:],color='black')
    cbar = plt.colorbar(cp)
    cbar.set_label(cbarLabel[k])
    if(showDensity and (dynName[k] == 'dynDiag')):
        salt = np.squeeze(np.average(data[1,:,1:-1,:],weights=openFrac[:,1:-1,:],axis=1))
        pressure = -1 * np.ones(salt.shape) * 1020 * 9.81 * np.repeat(np.expand_dims(z,1), salt.shape[1], axis=1) /10e3
        CT = gsw.CT_from_t(salt, np.average(data[0,:,1:-1,:],weights=openFrac[:,1:-1,:],axis=1), pressure)
        density = gsw.rho(salt, CT, 0) - 1000 #in-stu density less 1000
        densityLevels = np.linspace(25,28,16)
        cc = ax.contour(
            np.squeeze(x[0,:]),
            np.squeeze(z),
            np.squeeze(density),
            densityLevels,
            colors='black',
            linewidths=0.5,
            alpha=0.5
        )
        plt.clabel(cc, inline=3, fontsize=8)
    if(showZeros):
        if( k == 2 or k == 6 or k == 7):
            cc = ax.contour(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(plotData),
                # np.squeeze(data[kk, :, 7, :]),
                [0],
                colors='gray',
                linewidths=0.5,
                alpha=0.5
            )
            plt.clabel(cc, inline=3, fontsize=8)
    ax.set_title(f"{name[k]} at {int(i*dt/86400)} days")
    ## Berg shading
    if(showBergs):
        bergShading = dataBerg[2,:,:,:]
        bergShading[bergShading > 0] = 1
        cp2 = ax.contourf(
        x[0,:],
        np.squeeze(z),
        np.squeeze(np.average(bergShading,axis=1)),
        [.25,.5,.75,.9],
        extend="max",
        alpha=.2,
        cmap='cmo.gray_r')
    # cbar2 = plt.colorbar(cp2)
    # cbar2.set_label('Ocean Fraction')
    ax.set_xlim([0,50])        
    if(ax == ax2):
        ax.set_xlabel('Along Fjord [km]')
    else:
        pass
        # ax.set_xticks([])
    # ax.set_xlabel('Along Fjord [m] min: %.2f %s, max: %.2f %s' %(np.nanmin(data[kk, :, 1:-1, :]),cbarLabel[k],np.nanmax(data[kk, :, 1:-1, :]),cbarLabel[k]))
    ax.set_ylabel('Depth [m]')

## Now we plot the Melange over time

# print(args)
fileStr = args.files
files = sorted(glob.glob('%s*.pickle'%fileStr))

toIterate = np.arange(len(files))
H0Time = np.zeros(np.shape(toIterate))
VTime = np.zeros(np.shape(toIterate))
UcTime = np.zeros(np.shape(toIterate))
lengthTime = np.zeros(np.shape(toIterate))
iterationNumber = np.zeros(np.shape(toIterate))
pressTime = np.zeros(np.shape(toIterate))
BTime = np.zeros(np.shape(toIterate))
gTime = np.zeros(np.shape(toIterate))
spdTime = np.zeros(np.shape(toIterate))
timeTime = np.zeros(np.shape(toIterate))
for j in toIterate:
    file = files[j]
    with open(files[j], 'rb') as fileName:
        data = pickle.load(fileName)
        fileName.close()
    it = it = file.replace('./', '').replace('.pickle', '').replace(fileStr,'')
    H0Time[j]=data.H0
    X_ = np.concatenate(([data.X[0]], data.X_, [data.X[-1]]))
    H = np.concatenate(([data.H0], data.H, [data.HL]))
    W = np.concatenate(([data.W0], data.W, [data.WL]))
    VTime[j] = simpson(H*W, x=X_)*1e-9
    lengthTime[j]=data.X[-1]
    UcTime[j] = data.Uc
    BTime[j] = -1 * np.mean(data.B) /365.0 # we flip this for plotting purposes
    timeTime[j] = data.t*365.0 
    pressTime[j] = data.H0*data.pressure(data.H0)
    spdTime[j]=np.mean(data.U)
    gTime[j]=np.mean(data.gg)
    iterationNumber[j]=int(it)
timeTime = timeTime - timeTime[0] # shift to start at t=0
timeTime = (timeTime - timeTime[0])/5 #if pre-code fix for total time

#This lets you swap out the x axis for melt/calving if you want
drivingVariable = timeTime
drivingLabel = 'Time [days]'
xString = 'Time'

ax3.plot(drivingVariable,pressTime,'-',color='xkcd:sunflower')
ax3.scatter(drivingVariable[0],pressTime[0],s=50,marker='*',color='black')
# sca=ax3.scatter(BTime,pressTime,s=None,c=timeTime,cmap='cividis')
# cbar=plt.colorbar(sca)
# cbar.set_label('Time [days]')
ax3.set_ylabel('Buttressing (Force/Width) [$Nm^{-1}$]')
ax3.set_yscale('log')
# ax3.set_xlim([-20,960])
ax3.set_xlabel(drivingLabel)
ax3.grid(alpha=.25)
ax3.set_title('Buttressing vs Driving Temp')
ax3_2 = ax3.twinx()
ax3_2.plot(mitgcmTime/86400.0,mitgcmTemp,'-',color='xkcd:dark rose')
# ax4_2.scatter(drivingVariable[0],H0Time[0],s=50,marker='*',color='black')
ax3_2.set_ylabel('Driving Temp [C]',color='xkcd:dark rose')
# ax3_2.tick_params(axis='y',labelcolor='xkcd:dark rose')
# ax3_2.grid(alpha=.25,color='xkcd:dark rose')


ax4.plot(drivingVariable,lengthTime,'-',color='xkcd:blueberry')
ax4.scatter(drivingVariable[0],lengthTime[0],s=50,marker='*',color='black')
ax4.set_ylabel('Mélange length [m]',color='xkcd:blueberry')
ax4.tick_params(axis='y',labelcolor='xkcd:blueberry')
ax4.grid(alpha=.25,color='xkcd:blueberry')
# ax4.set_xlim([-20,960])
# sca=ax2.scatter(BTime,lengthTime,s=None,c=timeTime,cmap='cividis')
# cbar=plt.colorbar(sca)
# cbar.set_label('Iteration')
ax4_2 = ax4.twinx()
ax4_2.plot(drivingVariable,H0Time,'-',color='xkcd:mulberry')
ax4_2.scatter(drivingVariable[0],H0Time[0],s=50,marker='*',color='black')
ax4_2.set_ylabel('Mélange H0 [m]',color='xkcd:mulberry')
ax4_2.tick_params(axis='y',labelcolor='xkcd:mulberry')
# ax4_2.grid(alpha=.25,color='xkcd:mulberry')
ax4.set_xlabel(drivingLabel)
ax4.set_title('Mélange Size')

sOption = ''
if(args.noDensity > 0):
    sOption = sOption + '_nD'
if(args.noZero > 0):
    sOption = sOption + '_nZ'
if(args.noBergs > 0):
    sOption = sOption + '_nB'
strTmp = f"figs/prettyPlot_{i:010}{sOption}.png"
# plt.tight_layout()
plt.savefig(strTmp, format='png', dpi=plotDPI)
plt.show()
plt.close()
