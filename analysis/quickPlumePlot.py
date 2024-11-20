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
    if "plumeDiag.0" in file:
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

os.system('rm -f figs/plumePlot*.png')
os.system('rm -f figs/plumePlot*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")



dynName = ['plumeDiag', 'plumeDiag', 'plumeDiag', 'plumeDiag','plumeDiag']
name = ["W", "Temp", "Sal", "Radius", "MeltRate"]
units = ["[m/s]", "[C]", "[PSU]", "[m]", "[m/day]"]


if(os.path.isfile('input/plumeMask.bin')):
    plumeMask = np.fromfile('input/plumeMask.bin', dtype='>f8')
    plumeMask = plumeMask.reshape(np.shape(x))
    plumeLocations = np.where(plumeMask == 3)

else:
    print('No Plumes Found, exiting PlumePlot')
    exit()

plumeLoc = [plumeLocations[0][0],plumeLocations[1][0]]
print('Plume Locaion is grid', plumeLoc )

for k in range(len(name)):
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        if k == 0:
            lvl = [1,3]
            cm = "xkcd:raspberry"
        elif k == 1:
            lvl = [33.5,35]
            cm = "xkcd:green"
        elif k == 2:
            lvl = [-0.5, 0.5]
            cm = "xkcd:rose"
        elif k == 3:
            lvl = [-0.5, 0.5]
            cm = "xkcd:violet"
        elif k == 4:
            lvl = [-0.005, 0.005]
            cm = "xkcd:lavender"
        elif k == 5:
            lvl = [-0.005, 0.005]
            cm = "xkcd:lavender"
        plt.figure()
        plt.plot(data[k,:,plumeLoc[0],plumeLoc[1]],np.squeeze(z),linewidth=1,color=cm)
        ax = plt.gca()
        # ax.set_xlim(lvl)

        plt.xlabel(name[k] + " " + units[k] + ' %.3f %.3f nan: %i' %(np.nanmin(data[k, :, plumeLoc[0],plumeLoc[1]]),np.nanmax(data[k, :, plumeLoc[0],plumeLoc[1]]),np.max(np.isnan(data[k, :, plumeLoc[0],plumeLoc[1]]))))
        plt.ylabel('Depth [m]')
        plt.title("%s x = %i at %.02f days" % (name[k], x[0,plumeLoc[1]], i/86400.0*dt))
        j = i/sizeStep + startStep
        
        str = "figs/plumePlot%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png',dpi=plotDPI)
        plt.close()
        plt.show()

    os.system('magick -delay %f figs/plumePlot%s*.png -colors 256 -depth 256 figs/plumePlot%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
if(cleanPNGs):
    os.system('rm -f figs/plumePlot*.png')
