from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import fileinput

# Plot settings from local helper
zDepth = -50
plotDPI = 100
cleanPNGs = True

directories = ['fjord_b0','fjord_b1','fjord_m0']

if(os.path.isfile(directories[0] + '/input/plotHelperLocal.py')):
    sys.path.append(directories[0] + '/input')
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
for line in fileinput.input(directories[0] + '/input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

#Find time steps to take
maxStep = 0
sizeStep = 1e10
startStep = 1e10

for file in os.listdir(directories[0] + '/results'):
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

# os.system('rm -f figs/plumePlot*.png')
# os.system('rm -f figs/plumePlot*.gif')

x = mds.rdmds(directories[0] + "/results/XC")
y = mds.rdmds(directories[0] + "/results/YC")
z = mds.rdmds(directories[0] + "/results/RC")



dynName = ['plumeDiag', 'plumeDiag', 'plumeDiag', 'plumeDiag','plumeDiag']
name = ["W", "Temp", "Sal", "Total Melt Rate", "Plume Melt Rate"]
units = ["[m/s]", "[C]", "[PSU]", "[m/day]", "[m/day]"]


if(os.path.isfile(directories[0] + '/input/plumeMask.bin')):
    plumeMask = np.fromfile(directories[0] + '/input/plumeMask.bin', dtype='>f8')
    plumeMask = plumeMask.reshape(np.shape(x))
    plumeLocations = np.where(plumeMask == 3)

else:
    print('No Plumes Found, exiting PlumePlot')
    exit()

plumeLoc = [plumeLocations[0][0],plumeLocations[1][0]]
print('Plume Locaion is grid', plumeLoc )

for k in range(len(name)):
    plt.figure()
    for folder in directories:
        i = maxStep
        data = mds.rdmds(folder + "/results/%s"%(dynName[k]), i)
        # data[:,:,plumeLoc[0],plumeLoc[1]] = np.nan
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
        
        totalMean = np.nanmean(data[k,:,:,plumeLoc[1]])
        pl = plt.plot(np.nanmean(data[k,:,:,plumeLoc[1]],axis=1),np.squeeze(z),linewidth=1,label=folder)
        plt.plot([totalMean,totalMean],[np.min(z),np.max(z)],linewidth=1,linestyle='--',color=pl[0].get_color())
        ax = plt.gca()
        # ax.set_xlim(lvl)

        plt.xlabel(name[k] + " " + units[k] + ' %.3f %.3f nan: %i' %(np.nanmin(data[k, :, plumeLoc[0],plumeLoc[1]]),np.nanmax(data[k, :, plumeLoc[0],plumeLoc[1]]),np.max(np.isnan(data[k, :, plumeLoc[0],plumeLoc[1]]))))
        plt.ylabel('Depth [m]')
        plt.title("Width Averaged %s at %.02f days" % (name[k], i/86400.0*dt))
        j = i/sizeStep + startStep
            
    str = "figs/iceFront%s%05i.png" % (name[k],j)
    plt.legend()    
    plt.savefig(str, format='png',dpi=plotDPI)
    plt.show()
    plt.close()
    


