from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput

# Plot settings from local helper
zDepth = -50
plotDPI = 100
cleanPNGs = True

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

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

#Find time steps to take
maxStep = 0
sizeStep = 1e10
startStep = 1e10

for file in os.listdir('results'):
    # print(file)
    if "BRGFlx.0" in file:
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

os.system('rm -f figs/bergMap*.png')
os.system('rm -f figs/bergMap*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")



dynName = ['BRGFlx', 'BRGFlx', 'BRGFlx', 'BRGFlx','BRGFlx']
name = ['BRGfwFlx','BRGhtFlx','BRGmltRt','BRG_TauX','BRG_TauY']
units = ["[kg/m^2/s]", "[W/m^2]", "[m/d]", "[kN/m^2]","[kN/m^2]"]

zSlice = np.argmin(np.abs(z[:,0,0]- zDepth))
print('depth is z =', z[zSlice,0,0], 'index', zSlice)

for k in range(len(name)):
    print('\t',name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        if(data.shape[0] < (k + 1) ):
            if(i == startStep):
                print("\t\t%s not available, skipping" % name[k])
            break
        if k == 0:
            lvl = np.linspace(0,0.1,64)
            cm = "cmo.deep"
        elif k == 1:
            lvl = np.linspace(0,1000,64)
            cm = "cmo.amp"
        elif k == 2:
            lvl = np.linspace(0,2,32)
            cm = "cmo.speed"
        elif k == 3:
            lvl = np.linspace(-10,10,127)
            data = data/1000 # Pa to kPa
            cm = "cmo.balance"
        elif k == 4:
            lvl = np.linspace(-10,10,127)
            data = data/1000 # Pa to kPa
            cm = "cmo.balance"
        plt.figure()
        cp = plt.contourf(
            np.squeeze(x),
            np.squeeze(y),
            np.squeeze(data[k, zSlice, :, :]),
            lvl,
            extend="both",
            cmap=cm,
        )
        cbar = plt.colorbar(cp)
        cbar.set_label(units[k])
        plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(data[k, zSlice, :, :]),np.nanmax(data[k, zSlice, :, :]),np.max(np.isnan(data[k, zSlice, :, :]))))
        plt.ylabel('Depth [m]')
        plt.title("%s depth %f at %.02f days" % (name[k], z[zSlice,0,0] ,i/86400.0*dt))
        j = i/sizeStep + startStep
        
        str = "figs/bergMap%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png',dpi=plotDPI)
        plt.close()
        plt.show()
    if(not data.shape[0] < (k + 1) ):    
        os.system('magick -delay %f figs/bergMap%s*.png -colors 256 -depth 256 figs/bergMap%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/bergMap*.png')
