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
print('Plot DPI:',plotDPI,'; clean PNGs:',cleanPNGs, '; usePcolor:', usePcolor)

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
    if "shelfDiag.0" in file:
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

os.system('rm -f figs/shelfMap*.png')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")



dynName = ['shelfDiag', 'shelfDiag', 'shelfDiag', 'shelfDiag']
name = ['SHIfwFlx','SHIhtFlx','SHI_TauX','SHI_TauY']
units = ["[m^3/s]", "[W/m^2]", "[N/m^2]","[N/m^2]"]

zSlice = np.argmin(np.abs(z[:,0,0]- zDepth))
print('Plots values across bottom of shelf')

for k in range(len(name)):
    print('\t',name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        if k == 0: #
            lvl = np.linspace(0,3,128)
            cm = "cmo.tempo"
            # raw data is kg/m^2/s, convert to m^3/s
            # kg / (m^2 s) * (m^3/1000 kg) * dx [m] * dy [m]
            plotData = -1 * data[k, :, :] / 1000 * (x[1,1] - x[1,0]) * (y[1,1] - y[0,1])
            # print((x[1,1] - x[1,0]),(y[1,1] - y[0,1]))
            # print(data.shape)
            # print(plotData.shape)
        elif k == 1:
            lvl = np.linspace(0,5,64)
            cm = "cmo.amp"
            plotData = data[k,:, :]
        elif k == 2:
            lvl = np.linspace(-1.5,1.5,127)
            plotData = data[k,:, :] # Pa 
            cm = "cmo.balance"
        elif k == 3:
            lvl = np.linspace(-1.5,1.5,127)
            plotData = data[k,:, :] # Pa 
            cm = "cmo.balance"
        plt.figure()
        if(usePcolor):
            cp = plt.pcolormesh(
                np.squeeze(x),
                np.squeeze(y),
                np.squeeze(plotData),
                cmap=cm,
                vmin=np.min(lvl),
                vmax=np.max(lvl),
            )
        else:
            cp = plt.contourf(
                np.squeeze(x),
                np.squeeze(y),
                np.squeeze(plotData),
                lvl,
                extend="both",
                cmap=cm,
            )
        cbar = plt.colorbar(cp)
        cbar.set_label(units[k])
        plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(plotData),np.nanmax(plotData),np.max(np.isnan(plotData))))
        plt.ylabel('Depth [m]')
        plt.title("%s across shelf at %.02f days" % (name[k], i/86400.0*dt))
        j = i/sizeStep + startStep
        
        str = "figs/shelfMap%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png',dpi=plotDPI)
        # plt.show()
        plt.close()  
    os.system('magick -delay %f figs/shelfMap%s*.png -colors 256 -depth 256 figs/shelfMap%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/shelfMap*.png')
